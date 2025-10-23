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
        const sfcgalClass = new Module.SFCGAL();
        sfcgalClass.initialize();

        sfcgalInstance = sfcgalClass;
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
 * Wrap a SFCGAL function call with error handling
 * Catches Wasm runtime errors and converts them to SFCGALError
 *
 * @param {Function} fn - Function to wrap
 * @param {string} operationName - Name of operation for error messages
 * @returns {Function} Wrapped function
 */
export function wrapSFCGALOperation(fn, operationName = 'operation') {
    return (...args) => {
        try {
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
                if (result === 'GEOMETRYCOLLECTION EMPTY') {
                    console.warn(`[SFCGAL] ${operationName} returned empty geometry - geometries may not intersect`);
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
