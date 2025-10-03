// Helper functions for 2D Canvas rendering of geometries

function parseWKT2D(wkt) {
    const coords = [];
    const matches = wkt.matchAll(/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/g);
    for (const m of matches) {
        coords.push({ x: parseFloat(m[1]), y: parseFloat(m[2]) });
    }
    return coords;
}

function getBounds(coords) {
    if (coords.length === 0) return { minX: 0, maxX: 10, minY: 0, maxY: 10 };
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    coords.forEach(c => {
        minX = Math.min(minX, c.x);
        maxX = Math.max(maxX, c.x);
        minY = Math.min(minY, c.y);
        maxY = Math.max(maxY, c.y);
    });
    const margin = Math.max(maxX - minX, maxY - minY) * 0.1 || 1;
    return { minX: minX - margin, maxX: maxX + margin, minY: minY - margin, maxY: maxY + margin };
}

function createCanvas2DRenderer(canvasId) {
    const canvas = document.getElementById(canvasId);
    const ctx = canvas.getContext('2d');
    
    // Set canvas size
    canvas.width = canvas.parentElement.clientWidth;
    canvas.height = canvas.parentElement.clientHeight;
    
    let bounds = { minX: -1, maxX: 11, minY: -1, maxY: 11 };
    
    function toScreen(x, y) {
        const w = canvas.width;
        const h = canvas.height;
        const rangeX = bounds.maxX - bounds.minX;
        const rangeY = bounds.maxY - bounds.minY;
        return {
            x: (x - bounds.minX) / rangeX * w,
            y: h - (y - bounds.minY) / rangeY * h  // Flip Y axis
        };
    }
    
    function clear() {
        ctx.fillStyle = '#1a1a2e';
        ctx.fillRect(0, 0, canvas.width, canvas.height);
    }
    
    function drawGrid() {
        ctx.strokeStyle = '#333355';
        ctx.lineWidth = 1;
        
        const rangeX = bounds.maxX - bounds.minX;
        const rangeY = bounds.maxY - bounds.minY;
        const stepX = Math.pow(10, Math.floor(Math.log10(rangeX))) || 1;
        const stepY = Math.pow(10, Math.floor(Math.log10(rangeY))) || 1;
        
        // Vertical lines
        for (let x = Math.floor(bounds.minX / stepX) * stepX; x <= bounds.maxX; x += stepX) {
            const p = toScreen(x, 0);
            ctx.beginPath();
            ctx.moveTo(p.x, 0);
            ctx.lineTo(p.x, canvas.height);
            ctx.stroke();
        }
        
        // Horizontal lines
        for (let y = Math.floor(bounds.minY / stepY) * stepY; y <= bounds.maxY; y += stepY) {
            const p = toScreen(0, y);
            ctx.beginPath();
            ctx.moveTo(0, p.y);
            ctx.lineTo(canvas.width, p.y);
            ctx.stroke();
        }
        
        // Axes
        ctx.strokeStyle = '#444466';
        ctx.lineWidth = 2;
        const origin = toScreen(0, 0);
        ctx.beginPath();
        ctx.moveTo(0, origin.y);
        ctx.lineTo(canvas.width, origin.y);
        ctx.stroke();
        ctx.beginPath();
        ctx.moveTo(origin.x, 0);
        ctx.lineTo(origin.x, canvas.height);
        ctx.stroke();
    }
    
    function drawPolygon(wkt, color, opacity = 1, filled = true) {
        const coords = parseWKT2D(wkt);
        if (coords.length < 2) return;
        
        ctx.globalAlpha = opacity;
        ctx.beginPath();
        const start = toScreen(coords[0].x, coords[0].y);
        ctx.moveTo(start.x, start.y);
        for (let i = 1; i < coords.length; i++) {
            const p = toScreen(coords[i].x, coords[i].y);
            ctx.lineTo(p.x, p.y);
        }
        ctx.closePath();
        
        if (filled) {
            ctx.fillStyle = color;
            ctx.fill();
        }
        
        ctx.strokeStyle = color;
        ctx.lineWidth = 2;
        ctx.stroke();
        ctx.globalAlpha = 1;
    }
    
    function drawPoint(x, y, color, radius = 5) {
        const p = toScreen(x, y);
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.arc(p.x, p.y, radius, 0, Math.PI * 2);
        ctx.fill();
    }
    
    function drawLine(x1, y1, x2, y2, color, lineWidth = 2, dashed = false) {
        const p1 = toScreen(x1, y1);
        const p2 = toScreen(x2, y2);
        ctx.strokeStyle = color;
        ctx.lineWidth = lineWidth;
        if (dashed) {
            ctx.setLineDash([5, 5]);
        } else {
            ctx.setLineDash([]);
        }
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.stroke();
        ctx.setLineDash([]);
    }
    
    function updateBounds(...wkts) {
        const allCoords = [];
        wkts.forEach(wkt => {
            if (wkt) allCoords.push(...parseWKT2D(wkt));
        });
        if (allCoords.length > 0) {
            bounds = getBounds(allCoords);
        }
    }
    
    function render() {
        clear();
        drawGrid();
    }
    
    return {
        ctx,
        canvas,
        clear,
        drawGrid,
        drawPolygon,
        drawPoint,
        drawLine,
        updateBounds,
        render,
        toScreen
    };
}

export { createCanvas2DRenderer, parseWKT2D, getBounds };
