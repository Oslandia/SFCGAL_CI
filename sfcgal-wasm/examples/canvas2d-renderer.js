/**
 * Canvas 2D Renderer for SFCGAL Geometries
 * Renders WKT geometries on HTML5 Canvas with proper scaling and colors
 */

export class Canvas2DRenderer {
    constructor(canvasId, options = {}) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) {
            throw new Error(`Canvas element #${canvasId} not found`);
        }

        this.ctx = this.canvas.getContext('2d');
        this.options = {
            padding: options.padding || 40,
            fillColor: options.fillColor || 'rgba(102, 126, 234, 0.3)',
            strokeColor: options.strokeColor || 'rgb(102, 126, 234)',
            lineWidth: options.lineWidth || 2,
            pointRadius: options.pointRadius || 4,
            gridColor: options.gridColor || '#e0e0e0',
            axisColor: options.axisColor || '#999',
            showGrid: options.showGrid !== false,
            showAxes: options.showAxes !== false,
            backgroundColor: options.backgroundColor || '#ffffff'
        };

        this.geometries = [];
        this.bounds = null;
    }

    clear() {
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        this.ctx.fillStyle = this.options.backgroundColor;
        this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
    }

    clearGeometries() {
        this.geometries = [];
        this.bounds = null;
    }

    addGeometry(wkt, options = {}) {
        this.geometries.push({ wkt, options });
    }

    render() {
        if (this.geometries.length === 0) return;

        // Calculate bounds
        this.calculateBounds();

        // Clear canvas
        this.clear();

        // Draw grid and axes
        if (this.options.showGrid) this.drawGrid();
        if (this.options.showAxes) this.drawAxes();

        // Draw each geometry
        for (const { wkt, options } of this.geometries) {
            this.drawWKT(wkt, options);
        }
    }

    calculateBounds() {
        let minX = Infinity, minY = Infinity;
        let maxX = -Infinity, maxY = -Infinity;

        for (const { wkt } of this.geometries) {
            const coords = this.extractCoords(wkt);
            for (const [x, y] of coords) {
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
            }
        }

        // Add margin
        const marginX = (maxX - minX) * 0.1 || 1;
        const marginY = (maxY - minY) * 0.1 || 1;

        this.bounds = {
            minX: minX - marginX,
            minY: minY - marginY,
            maxX: maxX + marginX,
            maxY: maxY + marginY
        };
    }

    worldToCanvas(x, y) {
        if (!this.bounds) return [0, 0];

        const { minX, minY, maxX, maxY } = this.bounds;
        const { padding } = this.options;

        const canvasWidth = this.canvas.width - 2 * padding;
        const canvasHeight = this.canvas.height - 2 * padding;

        const scaleX = canvasWidth / (maxX - minX);
        const scaleY = canvasHeight / (maxY - minY);
        const scale = Math.min(scaleX, scaleY);

        const offsetX = (canvasWidth - scale * (maxX - minX)) / 2;
        const offsetY = (canvasHeight - scale * (maxY - minY)) / 2;

        const cx = padding + offsetX + (x - minX) * scale;
        const cy = this.canvas.height - (padding + offsetY + (y - minY) * scale);

        return [cx, cy];
    }

    extractCoords(wkt) {
        const coords = [];
        const regex = /-?\d+\.?\d*\s+-?\d+\.?\d*/g;
        const matches = wkt.match(regex);

        if (matches) {
            for (const match of matches) {
                const [x, y] = match.trim().split(/\s+/).map(Number);
                coords.push([x, y]);
            }
        }

        return coords;
    }

    drawWKT(wkt, options = {}) {
        const opts = { ...this.options, ...options };

        this.ctx.strokeStyle = opts.strokeColor;
        this.ctx.fillStyle = opts.fillColor;
        this.ctx.lineWidth = opts.lineWidth;

        const upperWKT = wkt.toUpperCase();

        if (upperWKT.includes('POINT')) {
            this.drawPoint(wkt, opts);
        } else if (upperWKT.includes('LINESTRING') || upperWKT.includes('MULTILINESTRING')) {
            this.drawLineString(wkt, opts);
        } else if (upperWKT.includes('POLYGON') || upperWKT.includes('MULTIPOLYGON')) {
            this.drawPolygon(wkt, opts);
        } else if (upperWKT.includes('GEOMETRYCOLLECTION')) {
            // Extract sub-geometries (simplified)
            console.warn('GEOMETRYCOLLECTION rendering not fully implemented');
        }
    }

    drawPoint(wkt, opts) {
        const coords = this.extractCoords(wkt);
        this.ctx.fillStyle = opts.strokeColor;

        for (const [x, y] of coords) {
            const [cx, cy] = this.worldToCanvas(x, y);
            this.ctx.beginPath();
            this.ctx.arc(cx, cy, opts.pointRadius, 0, 2 * Math.PI);
            this.ctx.fill();
        }
    }

    drawLineString(wkt, opts) {
        const coords = this.extractCoords(wkt);
        if (coords.length < 2) return;

        this.ctx.beginPath();
        const [startX, startY] = this.worldToCanvas(coords[0][0], coords[0][1]);
        this.ctx.moveTo(startX, startY);

        for (let i = 1; i < coords.length; i++) {
            const [cx, cy] = this.worldToCanvas(coords[i][0], coords[i][1]);
            this.ctx.lineTo(cx, cy);
        }

        this.ctx.stroke();
    }

    drawPolygon(wkt, opts) {
        // Extract outer ring (simplified - assumes single polygon)
        const coords = this.extractCoords(wkt);
        if (coords.length < 3) return;

        this.ctx.beginPath();
        const [startX, startY] = this.worldToCanvas(coords[0][0], coords[0][1]);
        this.ctx.moveTo(startX, startY);

        for (let i = 1; i < coords.length; i++) {
            const [cx, cy] = this.worldToCanvas(coords[i][0], coords[i][1]);
            this.ctx.lineTo(cx, cy);
        }

        this.ctx.closePath();
        this.ctx.fill();
        this.ctx.stroke();
    }

    drawGrid() {
        if (!this.bounds) return;

        this.ctx.strokeStyle = this.options.gridColor;
        this.ctx.lineWidth = 1;

        const { minX, minY, maxX, maxY } = this.bounds;
        const steps = 10;
        const stepX = (maxX - minX) / steps;
        const stepY = (maxY - minY) / steps;

        // Vertical lines
        for (let i = 0; i <= steps; i++) {
            const x = minX + i * stepX;
            const [cx1, cy1] = this.worldToCanvas(x, minY);
            const [cx2, cy2] = this.worldToCanvas(x, maxY);
            this.ctx.beginPath();
            this.ctx.moveTo(cx1, cy1);
            this.ctx.lineTo(cx2, cy2);
            this.ctx.stroke();
        }

        // Horizontal lines
        for (let i = 0; i <= steps; i++) {
            const y = minY + i * stepY;
            const [cx1, cy1] = this.worldToCanvas(minX, y);
            const [cx2, cy2] = this.worldToCanvas(maxX, y);
            this.ctx.beginPath();
            this.ctx.moveTo(cx1, cy1);
            this.ctx.lineTo(cx2, cy2);
            this.ctx.stroke();
        }
    }

    drawAxes() {
        if (!this.bounds) return;

        this.ctx.strokeStyle = this.options.axisColor;
        this.ctx.lineWidth = 2;

        const { minX, minY, maxX, maxY } = this.bounds;

        // X axis
        const [x1, y1] = this.worldToCanvas(minX, 0);
        const [x2, y2] = this.worldToCanvas(maxX, 0);
        this.ctx.beginPath();
        this.ctx.moveTo(x1, y1);
        this.ctx.lineTo(x2, y2);
        this.ctx.stroke();

        // Y axis
        const [x3, y3] = this.worldToCanvas(0, minY);
        const [x4, y4] = this.worldToCanvas(0, maxY);
        this.ctx.beginPath();
        this.ctx.moveTo(x3, y3);
        this.ctx.lineTo(x4, y4);
        this.ctx.stroke();
    }
}
