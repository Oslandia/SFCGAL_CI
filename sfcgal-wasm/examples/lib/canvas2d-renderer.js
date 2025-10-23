/**
 * Canvas 2D Renderer for SFCGAL geometries
 * Renders WKT geometries on HTML5 Canvas with automatic scaling and grid display
 *
 * @class Canvas2DRenderer
 * @example
 * const renderer = new Canvas2DRenderer('myCanvas', { showGrid: true });
 * await renderer.addGeometry('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))', { color: '#3b82f6' });
 */
export class Canvas2DRenderer {
    /**
     * Create a new Canvas2D renderer
     * @param {string} canvasId - ID of the canvas element
     * @param {Object} [options={}] - Rendering options
     * @param {boolean} [options.showGrid=true] - Show coordinate grid
     * @param {boolean} [options.showAxes=true] - Show X/Y axes
     * @throws {Error} If canvas element is not found
     */
    constructor(canvasId, options = {}) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) throw new Error(`Canvas '${canvasId}' not found`);

        this.ctx = this.canvas.getContext('2d');
        this.geometries = [];
        this.colors = ['#ef4444', '#3b82f6', '#10b981', '#f59e0b', '#8b5cf6'];
        this.padding = 20;
        this.options = { showGrid: true, showAxes: true, ...options };

        // Setup high DPI canvas
        this.setupHighDPI();
    }

    /**
     * Setup high DPI (Retina) canvas rendering
     * Adjusts canvas dimensions based on device pixel ratio for sharp rendering
     * @private
     */
    setupHighDPI() {
        const dpr = window.devicePixelRatio || 1;
        const rect = this.canvas.getBoundingClientRect();

        this.canvas.width = rect.width * dpr;
        this.canvas.height = rect.height * dpr;
        this.canvas.style.width = rect.width + 'px';
        this.canvas.style.height = rect.height + 'px';

        this.ctx.scale(dpr, dpr);
        this.displayWidth = rect.width;
        this.displayHeight = rect.height;
    }

    /**
     * Parse WKT (Well-Known Text) string into internal geometry representation
     * @param {string} wkt - WKT geometry string
     * @returns {Object} Parsed geometry object with type and coordinates
     * @returns {string} return.type - Geometry type (POINT, LINESTRING, POLYGON, etc.)
     * @returns {Array} return.coords - Array of coordinate pairs [[x, y], ...]
     * @throws {Error} If WKT is empty, invalid, or unsupported
     * @private
     */
    parseWKT(wkt) {
        wkt = wkt.trim();

        // Handle empty or invalid WKT
        if (!wkt || wkt === 'GEOMETRYCOLLECTION EMPTY' || wkt.includes('EMPTY')) {
            throw new Error(`Empty or invalid WKT: ${wkt}`);
        }

        // POINT
        if (wkt.match(/^POINT\s*\(/i)) {
            const match = wkt.match(/POINT\s*\(([^)]+)\)/i);
            if (match) {
                const parts = match[1].trim().split(/\s+/);
                const [x, y] = parts.map(Number);
                if (isNaN(x) || isNaN(y)) {
                    throw new Error(`Invalid POINT coordinates: ${match[1]}`);
                }
                return { type: 'POINT', coords: [[x, y]] };
            }
        }

        // LINESTRING
        if (wkt.match(/^LINESTRING\s*\(/i)) {
            const match = wkt.match(/LINESTRING\s*\(([^)]+)\)/i);
            if (match) {
                const coords = match[1].split(',').map(pair => {
                    const [x, y] = pair.trim().split(/\s+/).map(Number);
                    return [x, y];
                }).filter(([x, y]) => !isNaN(x) && !isNaN(y));

                if (coords.length === 0) {
                    throw new Error('Invalid LINESTRING: no valid coordinates found');
                }
                return { type: 'LINESTRING', coords };
            }
        }

        // TRIANGLE (same as POLYGON but with specific name)
        if (wkt.match(/^TRIANGLE\s*\(/i)) {
            const match = wkt.match(/TRIANGLE\s*\(\(([^)]+)\)\)/i);
            if (match) {
                const coords = match[1].split(',').map(pair => {
                    const [x, y] = pair.trim().split(/\s+/).map(Number);
                    return [x, y];
                }).filter(([x, y]) => !isNaN(x) && !isNaN(y));

                if (coords.length < 3) {
                    throw new Error(`Invalid TRIANGLE: requires at least 3 points, got ${coords.length}`);
                }
                return { type: 'POLYGON', coords };
            }
        }

        // POLYGON - improved to handle more cases
        if (wkt.match(/^POLYGON\s*\(/i)) {
            // Try to match simple polygon: POLYGON((coords))
            let match = wkt.match(/POLYGON\s*\(\s*\(([^)]+)\)\s*\)/i);
            if (match) {
                const coords = match[1].split(',').map(pair => {
                    const [x, y] = pair.trim().split(/\s+/).map(Number);
                    return [x, y];
                }).filter(([x, y]) => !isNaN(x) && !isNaN(y));

                if (coords.length < 3) {
                    throw new Error(`Invalid POLYGON: requires at least 3 points, got ${coords.length}`);
                }
                return { type: 'POLYGON', coords };
            }
        }

        // MULTIPOINT - Robust parser for both formats
        if (wkt.match(/^MULTIPOINT\s*\(/i)) {
            // Parsing MULTIPOINT
            const coords = [];

            // Remove MULTIPOINT prefix and outer parentheses
            const content = wkt.replace(/^MULTIPOINT\s*\(/i, '').replace(/\)\s*$/, '');

            // Split by comma to get individual points
            const points = content.split(',');

            for (const point of points) {
                // Remove parentheses if present: (x y) -> x y
                const cleaned = point.trim().replace(/^\(/, '').replace(/\)$/, '');

                // Split by whitespace
                const parts = cleaned.split(/\s+/).filter(p => p.length > 0);

                if (parts.length >= 2) {
                    const x = parseFloat(parts[0]);
                    const y = parseFloat(parts[1]);

                    if (!isNaN(x) && !isNaN(y)) {
                        coords.push([x, y]);
                    } else {
                        // Skip invalid coordinates
                    }
                }
            }

            if (coords.length === 0) {
                throw new Error('No valid MULTIPOINT coordinates found');
            }

            // MULTIPOINT parsed successfully
            return { type: 'MULTIPOINT', coords };
        }

        // MULTILINESTRING
        if (wkt.match(/^MULTILINESTRING\s*\(/i)) {
            const lines = [];
            const matches = [...wkt.matchAll(/\(([^)]+)\)/g)];
            for (const match of matches) {
                const coords = match[1].split(',').map(pair => {
                    const [x, y] = pair.trim().split(/\s+/).map(Number);
                    return [x, y];
                });
                lines.push(coords);
            }
            return { type: 'MULTILINESTRING', lines };
        }

        // MULTIPOLYGON
        if (wkt.match(/^MULTIPOLYGON\s*\(/i)) {
            const polygons = [];
            const polygonMatches = [...wkt.matchAll(/\(\(([^)]+)\)\)/g)];
            for (const match of polygonMatches) {
                const coords = match[1].split(',').map(pair => {
                    const [x, y] = pair.trim().split(/\s+/).map(Number);
                    return [x, y];
                });
                polygons.push(coords);
            }
            return { type: 'MULTIPOLYGON', polygons };
        }

        // GEOMETRYCOLLECTION - try to extract simple geometries
        if (wkt.match(/^GEOMETRYCOLLECTION\s*\(/i)) {
            return { type: 'GEOMETRYCOLLECTION', wkt };
        }

        throw new Error(`Unsupported WKT: ${wkt.substring(0, 30)}`);
    }

    getBounds() {
        if (this.geometries.length === 0) {
            return { minX: 0, minY: 0, maxX: 100, maxY: 100 };
        }

        let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;

        const updateBounds = (coords) => {
            coords.forEach(([x, y]) => {
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
            });
        };

        this.geometries.forEach(({ geom }) => {
            // Skip empty geometries
            if (geom.type === 'EMPTY') return;

            if (geom.coords) {
                updateBounds(geom.coords);
            }
            if (geom.lines) {
                geom.lines.forEach(updateBounds);
            }
            if (geom.polygons) {
                geom.polygons.forEach(updateBounds);
            }
        });

        // Check if we found any valid bounds
        if (minX === Infinity || maxX === -Infinity) {
            return { minX: 0, minY: 0, maxX: 100, maxY: 100 };
        }

        const rangeX = maxX - minX || 10;
        const rangeY = maxY - minY || 10;
        const pad = 0.1;

        return {
            minX: minX - rangeX * pad,
            minY: minY - rangeY * pad,
            maxX: maxX + rangeX * pad,
            maxY: maxY + rangeY * pad
        };
    }

    worldToCanvas(x, y, bounds) {
        const width = this.displayWidth - 2 * this.padding;
        const height = this.displayHeight - 2 * this.padding;
        const worldWidth = bounds.maxX - bounds.minX;
        const worldHeight = bounds.maxY - bounds.minY;
        const scale = Math.min(width / worldWidth, height / worldHeight);
        const offsetX = (width - worldWidth * scale) / 2;
        const offsetY = (height - worldHeight * scale) / 2;

        const canvasX = this.padding + offsetX + (x - bounds.minX) * scale;
        const canvasY = this.displayHeight - (this.padding + offsetY + (y - bounds.minY) * scale);
        return [canvasX, canvasY];
    }

    drawGeometry(geom, opts, bounds) {
        // Skip empty geometries
        if (geom.type === 'EMPTY') {
            return;
        }

        const color = opts.color || opts.fillColor || '#3b82f6';
        const strokeColor = opts.strokeColor || color;
        const fillColor = opts.fillColor || this.hexToRGBA(color, opts.opacity || 0.3);
        const lineWidth = opts.lineWidth || 2;
        const pattern = opts.pattern || 'solid';

        this.ctx.strokeStyle = strokeColor;
        this.ctx.lineWidth = lineWidth;

        const drawPath = (coords) => {
            this.ctx.beginPath();
            coords.forEach(([x, y], i) => {
                const [cx, cy] = this.worldToCanvas(x, y, bounds);
                if (i === 0) this.ctx.moveTo(cx, cy);
                else this.ctx.lineTo(cx, cy);
            });
        };

        if (geom.type === 'POINT') {
            const [px, py] = this.worldToCanvas(geom.coords[0][0], geom.coords[0][1], bounds);
            this.ctx.fillStyle = color;
            this.ctx.beginPath();
            this.ctx.arc(px, py, opts.pointRadius || 6, 0, Math.PI * 2);
            this.ctx.fill();
            this.ctx.stroke();
        } else if (geom.type === 'LINESTRING') {
            drawPath(geom.coords);
            this.ctx.stroke();
        } else if (geom.type === 'POLYGON') {
            drawPath(geom.coords);
            this.ctx.closePath();

            // Apply fill pattern
            if (pattern === 'stripe' || pattern === 'hatch') {
                this.ctx.save();
                this.ctx.clip();
                this.drawStripePattern(bounds, strokeColor);
                this.ctx.restore();
            } else {
                this.ctx.fillStyle = fillColor;
                this.ctx.fill();
            }
            this.ctx.stroke();
        } else if (geom.type === 'MULTIPOINT') {
            this.ctx.fillStyle = color;
            geom.coords.forEach(([x, y]) => {
                const [px, py] = this.worldToCanvas(x, y, bounds);
                this.ctx.beginPath();
                this.ctx.arc(px, py, opts.pointRadius || 6, 0, Math.PI * 2);
                this.ctx.fill();
                this.ctx.stroke();
            });
        } else if (geom.type === 'MULTILINESTRING') {
            geom.lines.forEach(line => {
                drawPath(line);
                this.ctx.stroke();
            });
        } else if (geom.type === 'MULTIPOLYGON') {
            geom.polygons.forEach(poly => {
                drawPath(poly);
                this.ctx.closePath();
                if (pattern === 'stripe' || pattern === 'hatch') {
                    this.ctx.save();
                    this.ctx.clip();
                    this.drawStripePattern(bounds, strokeColor);
                    this.ctx.restore();
                } else {
                    this.ctx.fillStyle = fillColor;
                    this.ctx.fill();
                }
                this.ctx.stroke();
            });
        } else if (geom.type === 'GEOMETRYCOLLECTION') {
            this.ctx.fillStyle = '#999';
            this.ctx.font = '12px monospace';
            this.ctx.fillText('GEOMETRYCOLLECTION', 10, 20);
        }
    }

    hexToRGBA(hex, alpha) {
        const r = parseInt(hex.slice(1, 3), 16);
        const g = parseInt(hex.slice(3, 5), 16);
        const b = parseInt(hex.slice(5, 7), 16);
        return `rgba(${r}, ${g}, ${b}, ${alpha})`;
    }

    drawStripePattern(bounds, color) {
        this.ctx.strokeStyle = color;
        this.ctx.lineWidth = 1;
        const spacing = 8;
        for (let i = -this.displayHeight; i < this.displayWidth + this.displayHeight; i += spacing) {
            this.ctx.beginPath();
            this.ctx.moveTo(i, 0);
            this.ctx.lineTo(i + this.displayHeight, this.displayHeight);
            this.ctx.stroke();
        }
    }

    /**
     * Add a geometry to the canvas for rendering
     * @param {string} wkt - WKT geometry string
     * @param {Object|string} [optsOrColor=null] - Rendering options or color string
     * @param {string} [optsOrColor.color] - Fill/stroke color (hex)
     * @param {number} [optsOrColor.opacity=0.6] - Fill opacity (0-1)
     * @returns {Promise<void>}
     * @throws {Error} If WKT is invalid or unsupported
     * @example
     * // Using color string
     * renderer.addGeometry('POINT(5 5)', '#ff0000');
     *
     * // Using options object
     * renderer.addGeometry('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))', {
     *   color: '#3b82f6',
     *   opacity: 0.8
     * });
     */
    addGeometry(wkt, optsOrColor = null) {
        const geom = this.parseWKT(wkt);

        // Support both old API (color) and new API (options)
        let opts = {};
        if (typeof optsOrColor === 'string') {
            opts = { color: optsOrColor };
        } else if (optsOrColor && typeof optsOrColor === 'object') {
            opts = optsOrColor;
        }

        // Default color if not provided
        if (!opts.color && !opts.fillColor) {
            opts.color = this.colors[this.geometries.length % this.colors.length];
        }

        this.geometries.push({ geom, opts });
    }

    /**
     * Remove all geometries from the renderer
     * Does not redraw the canvas - call render() to update display
     * @returns {void}
     */
    clearGeometries() {
        this.geometries = [];
    }

    /**
     * Clear all geometries and redraw empty canvas
     * @returns {void}
     */
    clear() {
        this.ctx.clearRect(0, 0, this.displayWidth, this.displayHeight);

        if (this.options.showGrid) {
            this.ctx.strokeStyle = '#e5e7eb';
            this.ctx.lineWidth = 1;
            const gridSize = 40;
            for (let x = 0; x < this.displayWidth; x += gridSize) {
                this.ctx.beginPath();
                this.ctx.moveTo(x, 0);
                this.ctx.lineTo(x, this.displayHeight);
                this.ctx.stroke();
            }
            for (let y = 0; y < this.displayHeight; y += gridSize) {
                this.ctx.beginPath();
                this.ctx.moveTo(0, y);
                this.ctx.lineTo(this.displayWidth, y);
                this.ctx.stroke();
            }
        }
    }

    /**
     * Render all added geometries to the canvas
     * Automatically calculates bounds, applies scaling, and draws grid/axes if enabled
     * @returns {void}
     * @example
     * renderer.addGeometry('POINT(5 5)');
     * renderer.addGeometry('POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))');
     * renderer.render(); // Draw both geometries
     */
    render() {
        this.clear();
        if (this.geometries.length === 0) return;
        const bounds = this.getBounds();
        this.geometries.forEach(({ geom, opts }) => {
            this.drawGeometry(geom, opts, bounds);
        });
    }
}
