// Operation configurations and handlers

export const operations = {
    'create-geometry': {
        info: 'Create basic geometric shapes using WKT notation or predefined examples.',
        controls: `
            <div class="form-group">
                <label>Example Shapes</label>
                <div class="examples">
                    <button class="example-btn" onclick="loadExample('square')">Square</button>
                    <button class="example-btn" onclick="loadExample('hexagon')">Hexagon</button>
                    <button class="example-btn" onclick="loadExample('star')">Star</button>
                    <button class="example-btn" onclick="loadExample('circle')">Circle</button>
                </div>
            </div>

            <div class="form-group">
                <label>WKT Geometry</label>
                <textarea id="input-wkt" placeholder="POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <button class="btn" onclick="executeCreateGeometry()">Create & Display</button>

            <div id="result-create"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('input-wkt').value;
            if (!window.SFCGAL || !window.SFCGAL.isValid(wkt)) {
                alert('Invalid WKT geometry');
                return;
            }

            window.updateVisualization(wkt, 0);

            const analysis = window.SFCGAL.analyze(wkt);
            document.getElementById('result-create').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Geometry Information</div>
                    <div class="result-content">Type: ${analysis.type}
Dimension: ${analysis.dimension}
3D: ${analysis.is3D}
Valid: ${analysis.isValid}
Points: ${analysis.numPoints || 'N/A'}</div>
                </div>
            `;
        }
    },

    'extrude': {
        info: 'Extrude a 2D polygon into a 3D solid by adding height.',
        controls: `
            <div class="form-group">
                <label>Base Polygon (WKT)</label>
                <textarea id="extrude-input">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <div class="form-group">
                <label>Height: <span class="range-value" id="extrude-height-value">5</span></label>
                <input type="range" id="extrude-height" min="0.5" max="20" value="5" step="0.5"
                       oninput="document.getElementById('extrude-height-value').textContent = this.value">
            </div>

            <button class="btn" onclick="executeExtrude()">Extrude</button>

            <div id="result-extrude"></div>

            <div class="stats-grid" id="extrude-stats"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('extrude-input').value;
            const height = parseFloat(document.getElementById('extrude-height').value);

            if (!window.SFCGAL) return;

            const result = window.SFCGAL.extrudeDetailed(wkt, height, {});

            window.updateVisualization(wkt, height);

            document.getElementById('result-extrude').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Extruded Geometry</div>
                    <div class="result-content">${result.geometry.substring(0, 100)}...</div>
                </div>
            `;

            document.getElementById('extrude-stats').innerHTML = `
                <div class="stat-item">
                    <div class="stat-value">${result.baseArea?.toFixed(2) || 0}</div>
                    <div class="stat-label">Base Area</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${result.volume?.toFixed(2) || 0}</div>
                    <div class="stat-label">Volume</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${result.stats?.numFaces || 0}</div>
                    <div class="stat-label">Faces</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${result.stats?.surfaceArea?.toFixed(2) || 0}</div>
                    <div class="stat-label">Surface Area</div>
                </div>
            `;
        }
    },

    'area': {
        info: 'Calculate area (2D) or volume (3D) of geometries.',
        controls: `
            <div class="form-group">
                <label>Input Geometry (WKT)</label>
                <textarea id="area-input">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <button class="btn" onclick="executeArea()">Calculate</button>

            <div class="stats-grid" id="area-stats"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('area-input').value;
            if (!window.SFCGAL) return;

            const area = window.SFCGAL.area(wkt);
            const perimeter = window.SFCGAL.perimeter(wkt);
            const volume = window.SFCGAL.volume(wkt);
            const centroid = window.SFCGAL.centroid(wkt);

            window.updateVisualization(wkt, 0);

            document.getElementById('area-stats').innerHTML = `
                <div class="stat-item">
                    <div class="stat-value">${area.toFixed(2)}</div>
                    <div class="stat-label">Area</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${perimeter.toFixed(2)}</div>
                    <div class="stat-label">Perimeter</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${volume.toFixed(2)}</div>
                    <div class="stat-label">Volume</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">${centroid.substring(0, 20)}...</div>
                    <div class="stat-label">Centroid</div>
                </div>
            `;
        }
    },

    'intersection': {
        info: 'Compute the intersection of two geometries.',
        controls: `
            <div class="form-group">
                <label>Geometry A</label>
                <textarea id="int-geom-a">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <div class="form-group">
                <label>Geometry B</label>
                <textarea id="int-geom-b">POLYGON((5 5, 15 5, 15 15, 5 15, 5 5))</textarea>
            </div>

            <button class="btn" onclick="executeIntersection()">Compute Intersection</button>

            <div id="result-intersection"></div>
        `,
        execute: function() {
            const geomA = document.getElementById('int-geom-a').value;
            const geomB = document.getElementById('int-geom-b').value;

            if (!window.SFCGAL) return;

            const result = window.SFCGAL.intersection(geomA, geomB);

            window.updateVisualization(geomA, 0);

            document.getElementById('result-intersection').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Intersection Result</div>
                    <div class="result-content">${result}</div>
                </div>
            `;
        }
    },

    'union': {
        info: 'Compute the union of two geometries.',
        controls: `
            <div class="form-group">
                <label>Geometry A</label>
                <textarea id="union-geom-a">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <div class="form-group">
                <label>Geometry B</label>
                <textarea id="union-geom-b">POLYGON((5 5, 15 5, 15 15, 5 15, 5 5))</textarea>
            </div>

            <button class="btn" onclick="executeUnion()">Compute Union</button>

            <div id="result-union"></div>
        `,
        execute: function() {
            const geomA = document.getElementById('union-geom-a').value;
            const geomB = document.getElementById('union-geom-b').value;

            if (!window.SFCGAL) return;

            const result = window.SFCGAL.union(geomA, geomB);

            window.updateVisualization(geomA, 0);

            document.getElementById('result-union').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Union Result</div>
                    <div class="result-content">${result}</div>
                </div>
            `;
        }
    },

    'difference': {
        info: 'Compute the difference between two geometries (A - B).',
        controls: `
            <div class="form-group">
                <label>Geometry A</label>
                <textarea id="diff-geom-a">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <div class="form-group">
                <label>Geometry B</label>
                <textarea id="diff-geom-b">POLYGON((5 5, 15 5, 15 15, 5 15, 5 5))</textarea>
            </div>

            <button class="btn" onclick="executeDifference()">Compute Difference</button>

            <div id="result-difference"></div>
        `,
        execute: function() {
            const geomA = document.getElementById('diff-geom-a').value;
            const geomB = document.getElementById('diff-geom-b').value;

            if (!window.SFCGAL) return;

            const result = window.SFCGAL.difference(geomA, geomB);

            window.updateVisualization(geomA, 0);

            document.getElementById('result-difference').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Difference Result (A - B)</div>
                    <div class="result-content">${result}</div>
                </div>
            `;
        }
    },

    'translate': {
        info: 'Translate a geometry by a given offset in X, Y, Z directions.',
        controls: `
            <div class="form-group">
                <label>Input Geometry</label>
                <textarea id="translate-input">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <div class="form-group">
                <label>Offset X</label>
                <input type="number" id="translate-x" value="5" step="0.5">
            </div>

            <div class="form-group">
                <label>Offset Y</label>
                <input type="number" id="translate-y" value="5" step="0.5">
            </div>

            <div class="form-group">
                <label>Offset Z</label>
                <input type="number" id="translate-z" value="0" step="0.5">
            </div>

            <button class="btn" onclick="executeTranslate()">Translate</button>

            <div id="result-translate"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('translate-input').value;
            const dx = parseFloat(document.getElementById('translate-x').value);
            const dy = parseFloat(document.getElementById('translate-y').value);
            const dz = parseFloat(document.getElementById('translate-z').value);

            if (!window.SFCGAL) return;

            const result = window.SFCGAL.translate(wkt, dx, dy, dz);

            window.updateVisualization(result, 0);

            document.getElementById('result-translate').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Translated Geometry</div>
                    <div class="result-content">${result}</div>
                </div>
            `;
        }
    },

    'validate': {
        info: 'Validate a geometry and fix common issues.',
        controls: `
            <div class="form-group">
                <label>Input Geometry</label>
                <textarea id="validate-input">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <button class="btn" onclick="executeValidate()">Validate</button>

            <div id="result-validate"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('validate-input').value;

            if (!window.SFCGAL) return;

            const isValid = window.SFCGAL.isValid(wkt);
            const message = window.SFCGAL.validationMessage(wkt);

            let resultHTML = `
                <div class="result-area">
                    <div class="result-title">Validation Result</div>
                    <div class="result-content">Valid: ${isValid}
Message: ${message}</div>
                </div>
            `;

            if (!isValid) {
                const fixed = window.SFCGAL.makeValid(wkt);
                resultHTML += `
                    <div class="result-area" style="margin-top: 10px;">
                        <div class="result-title">Fixed Geometry</div>
                        <div class="result-content">${fixed}</div>
                    </div>
                `;
            }

            document.getElementById('result-validate').innerHTML = resultHTML;

            window.updateVisualization(wkt, 0);
        }
    },

    'analyze': {
        info: 'Analyze a geometry and get detailed information.',
        controls: `
            <div class="form-group">
                <label>Input Geometry</label>
                <textarea id="analyze-input">POLYGON((0 0, 10 0, 10 10, 0 10, 0 0))</textarea>
            </div>

            <button class="btn" onclick="executeAnalyze()">Analyze</button>

            <div id="result-analyze"></div>
        `,
        execute: function() {
            const wkt = document.getElementById('analyze-input').value;

            if (!window.SFCGAL) return;

            const analysis = window.SFCGAL.analyze(wkt);

            let info = `Type: ${analysis.type}
Dimension: ${analysis.dimension}
3D: ${analysis.is3D}
Valid: ${analysis.isValid}`;

            if (analysis.numPoints !== undefined) info += `\nPoints: ${analysis.numPoints}`;
            if (analysis.area !== undefined) info += `\nArea: ${analysis.area.toFixed(2)}`;
            if (analysis.perimeter !== undefined) info += `\nPerimeter: ${analysis.perimeter.toFixed(2)}`;
            if (analysis.numHoles !== undefined) info += `\nHoles: ${analysis.numHoles}`;

            document.getElementById('result-analyze').innerHTML = `
                <div class="result-area">
                    <div class="result-title">Analysis Result</div>
                    <div class="result-content">${info}</div>
                </div>
            `;

            window.updateVisualization(wkt, 0);
        }
    }
};