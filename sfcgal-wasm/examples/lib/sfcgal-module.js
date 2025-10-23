/**
 * SFCGAL WebAssembly Module Loader
 * Provides robust initialization with timeout and proper error handling
 */

let sfcgalInstance = null;

/**
 * Custom error class for SFCGAL operations
 */
export class SFCGALError extends Error {
    constructor(message, code, originalError = null) {
        super(message);
        this.name = 'SFCGALError';
        this.code = code;
        this.originalError = originalError;
    }
}

/**
 * Load and initialize the SFCGAL WebAssembly module with robust error handling
 * @param {Object} options - Configuration options
 * @param {number} options.timeout - Timeout in milliseconds (default: 30000)
 * @param {string} options.wasmPath - Custom path to SFCGAL wasm file
 * @returns {Promise<Object>} The SFCGAL instance with all methods
 * @throws {SFCGALError} If initialization fails or times out
 *
 * @example
 * const sfcgal = await loadSFCGAL({ timeout: 10000 });
 * const union = sfcgal.union("POINT(0 0)", "POINT(1 1)");
 */
export async function loadSFCGAL(options = {}) {
    if (sfcgalInstance) {
        return sfcgalInstance;
    }

    const { timeout = 30000, wasmPath = '../sfcgal.wasm' } = options;

    // Setup timeout abort controller
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), timeout);

    try {
        const SFCGALModule = await import('../sfcgal.js');

        // Initialize module with runtime ready promise
        const Module = await new Promise((resolve, reject) => {
            const moduleConfig = {
                locateFile: (path) => {
                    if (path.endsWith('.wasm')) {
                        return wasmPath;
                    }
                    return path;
                },
                onRuntimeInitialized: function() {
                    resolve(this);
                },
                onAbort: (err) => {
                    reject(new SFCGALError(
                        'WebAssembly module initialization aborted',
                        'WASM_ABORT',
                        err
                    ));
                },
            };

            SFCGALModule.default(moduleConfig).catch(reject);
        });

        clearTimeout(timeoutId);

        // Verify module is ready
        if (!Module || !Module.SFCGAL) {
            throw new SFCGALError(
                'SFCGAL class not found in module',
                'MODULE_INVALID'
            );
        }

        // Create SFCGAL instance
        const rawInstance = new Module.SFCGAL();
        rawInstance.initialize();

        // Wrap all SFCGAL methods with error handling
        const wrappedInstance = {};
        const methodsToWrap = [
            // Validation
            'isValid', 'getGeometryInfo',
            // Metrics
            'area', 'area3D', 'volume', 'length', 'length3D',
            'perimeter', 'perimeter3D', 'distance', 'distance3D',
            // Geometric operations
            'centroid', 'convexHull', 'buffer',
            // Boolean operations 2D
            'intersection', 'union', 'difference',
            // Boolean operations 3D
            'intersection3D', 'union3D', 'difference3D',
            // 3D operations
            'extrude', 'extrudeDetailed', 'toSolid',
            // Transformations
            'translate', 'rotate', 'scale',
            // Utility
            'version'
        ];

        methodsToWrap.forEach(method => {
            if (typeof rawInstance[method] === 'function') {
                wrappedInstance[method] = wrapSFCGALOperation(
                    rawInstance[method].bind(rawInstance),
                    method
                );
            } else {
                // Copy non-function properties as-is
                wrappedInstance[method] = rawInstance[method];
            }
        });

        sfcgalInstance = wrappedInstance;
        return sfcgalInstance;

    } catch (error) {
        clearTimeout(timeoutId);

        // Check if aborted due to timeout
        if (controller.signal.aborted) {
            throw new SFCGALError(
                `SFCGAL Wasm load timeout (>${timeout}ms)`,
                'LOAD_TIMEOUT'
            );
        }

        // Re-throw if already SFCGALError
        if (error instanceof SFCGALError) {
            throw error;
        }

        // Wrap other errors
        throw new SFCGALError(
            `Failed to load SFCGAL: ${error.message}`,
            'LOAD_FAILED',
            error
        );
    }
}

/**
 * Display a status message
 */
export function showStatus(elementId, message, type = 'info') {
    const element = document.getElementById(elementId);
    if (!element) return;

    element.textContent = message;
    element.className = 'status';

    if (type === 'success') {
        element.classList.add('status-success');
        setTimeout(() => { element.style.display = 'none'; }, 3000);
    } else if (type === 'error') {
        element.classList.add('status-error');
    } else {
        element.classList.add('status-info');
    }
}

export function getSFCGAL() {
    return sfcgalInstance;
}

/**
 * Validate WKT (Well-Known Text) format
 * @param {string} wkt - WKT string to validate
 * @returns {boolean} True if WKT appears valid
 */
export function isValidWKT(wkt) {
    if (!wkt || typeof wkt !== 'string') {
        return false;
    }

    // WKT pattern: GEOMETRY_TYPE( optional Z/M )( coordinates )
    const WKT_PATTERN = /^(POINT|LINESTRING|POLYGON|MULTIPOINT|MULTILINESTRING|MULTIPOLYGON|GEOMETRYCOLLECTION|POLYHEDRALSURFACE|SOLID|TRIANGULATEDSURFACE|TRIANGLE|TIN)\s*(\w+)?\s*\(/i;

    return WKT_PATTERN.test(wkt.trim());
}

/**
 * Validate WKT and throw error if invalid
 * @param {string} wkt - WKT string to validate
 * @param {string} operationName - Name of operation for error message
 * @throws {SFCGALError} If WKT is invalid
 */
export function validateWKT(wkt, operationName = 'operation') {
    if (!isValidWKT(wkt)) {
        throw new SFCGALError(
            `Invalid WKT format for ${operationName}. Expected format: GEOMETRY_TYPE(coordinates)`,
            'INVALID_WKT'
        );
    }
}

/**
 * Extract geometry type from WKT
 * @param {string} wkt - WKT string
 * @returns {string|null} Geometry type (e.g., 'POLYGON', 'SOLID') or null if invalid
 */
export function getGeometryType(wkt) {
    if (!wkt || typeof wkt !== 'string') return null;

    const match = wkt.trim().match(/^(\w+)\s*/i);
    return match ? match[1].toUpperCase() : null;
}

/**
 * Check if geometry is 3D (has Z coordinate)
 * @param {string} wkt - WKT string
 * @returns {boolean} True if geometry appears to be 3D
 */
export function is3DGeometry(wkt) {
    if (!wkt || typeof wkt !== 'string') return false;

    const type = getGeometryType(wkt);
    // Check for Z suffix or 3D geometry types
    return wkt.includes(' Z ') ||
           wkt.match(/\w+\s+Z\s*\(/i) !== null ||
           type === 'SOLID' ||
           type === 'POLYHEDRALSURFACE' ||
           type === 'TIN';
}

/**
 * Check if geometry type is suitable for volume calculation
 * @param {string} wkt - WKT string
 * @returns {boolean} True if geometry can have volume
 */
export function canCalculateVolume(wkt) {
    const type = getGeometryType(wkt);
    return type === 'SOLID' || type === 'POLYHEDRALSURFACE';
}

/**
 * Wrap a SFCGAL function call with error handling
 * Catches Wasm runtime errors and converts them to SFCGALError
 *
 * @param {Function} fn - Function to wrap
 * @param {string} operationName - Name of operation for error messages
 * @param {boolean} validateInput - Whether to validate WKT inputs (default: true)
 * @returns {Function} Wrapped function
 */
export function wrapSFCGALOperation(fn, operationName = 'operation', validateInput = true) {
    return (...args) => {
        try {
            // Validate WKT inputs if enabled
            if (validateInput) {
                args.forEach((arg, index) => {
                    if (typeof arg === 'string' && arg.trim().match(/^(POINT|LINESTRING|POLYGON|MULTI|GEOMETRYCOLLECTION|POLYHEDRALSURFACE|SOLID|TRIANGLE)/i)) {
                        validateWKT(arg, `${operationName} argument ${index + 1}`);
                    }
                });
            }

            const result = fn(...args);

            // Check for error strings returned by C++
            if (typeof result === 'string') {
                if (result.startsWith('ERROR:')) {
                    throw new SFCGALError(
                        result.substring(6),
                        'SFCGAL_ERROR'
                    );
                }

                // Check for EMPTY results that might indicate failure
                if (result === 'GEOMETRYCOLLECTION EMPTY' || result === 'POINT EMPTY' || result.includes(' EMPTY')) {
                    const emptyWarning = `${operationName} returned empty geometry - geometries may not intersect or may be degenerate`;
                    console.warn(`[SFCGAL] ${emptyWarning}`);

                    // Call global EMPTY callback if registered
                    if (typeof window !== 'undefined' && typeof window.onSFCGALEmptyResult === 'function') {
                        window.onSFCGALEmptyResult({
                            operation: operationName,
                            result: result,
                            message: emptyWarning
                        });
                    }
                }
            }

            return result;
        } catch (error) {
            // Catch Wasm runtime errors
            const errorMsg = error.message || String(error);

            if (errorMsg.includes('Aborted') ||
                errorMsg.includes('Cannot allocate memory') ||
                errorMsg.includes('Cannot enlarge memory') ||
                errorMsg.includes('divide by zero')) {
                throw new SFCGALError(
                    `${operationName} failed: Memory allocation error. The geometry may be too complex or invalid. Try simplifying the geometry or check its validity with isValid().`,
                    'WASM_ABORT',
                    error
                );
            }

            if (error instanceof SFCGALError) {
                throw error;
            }

            throw new SFCGALError(
                `${operationName} failed: ${errorMsg}`,
                'OPERATION_FAILED',
                error
            );
        }
    };
}
