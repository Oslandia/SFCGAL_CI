/**
 * Three.js 3D Renderer for SFCGAL geometries
 * Simplified version that actually works
 */

export class Three3DRenderer {
    constructor(containerId, options = {}) {
        this.containerId = containerId;
        this.container = document.getElementById(containerId);
        if (!this.container) {
            throw new Error(`Container '${containerId}' not found`);
        }

        this.THREE = null;
        this.OrbitControls = null;
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.geometries = [];
        this.colors = [0xef4444, 0x3b82f6, 0x10b981, 0xf59e0b, 0x8b5cf6];
        this.options = { wireframe: false, showGrid: true, showAxes: true, ...options };
        this.initialized = false;
        this.initPromise = this.initAsync();
    }

    async initAsync() {
        try {
            const [THREE, { OrbitControls }] = await Promise.all([
                import('three'),
                import('three/addons/controls/OrbitControls.js')
            ]);

            this.THREE = THREE;
            this.OrbitControls = OrbitControls;

            // Wait a bit for container to have dimensions
            await new Promise(resolve => {
                setTimeout(() => {
                    const width = this.container.clientWidth;
                    const height = this.container.clientHeight;

                    if (width > 0 && height > 0) {
                        this.init();
                        this.initialized = true;
                    } else {
                        console.warn('[Three3DRenderer] Container has zero dimensions');
                        this.initialized = false;
                    }
                    resolve();
                }, 100);
            });
        } catch (error) {
            console.error('Failed to load Three.js:', error);
            this.initialized = false;
        }
    }

    init() {
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;

        console.log('[Three3DRenderer] Init with dimensions:', width, 'x', height);

        // Scene
        this.scene = new this.THREE.Scene();
        this.scene.background = new this.THREE.Color(0xf0f0f0);

        // Camera
        this.camera = new this.THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        this.camera.position.set(10, 10, 10);
        this.camera.lookAt(0, 0, 0);

        // Renderer
        this.renderer = new this.THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(width, height);
        this.renderer.setPixelRatio(window.devicePixelRatio);

        // Clear and append
        this.container.innerHTML = '';
        this.container.appendChild(this.renderer.domElement);

        // Lights
        const ambientLight = new this.THREE.AmbientLight(0x404040, 2);
        this.scene.add(ambientLight);

        const directionalLight = new this.THREE.DirectionalLight(0xffffff, 1);
        directionalLight.position.set(10, 10, 10);
        this.scene.add(directionalLight);

        // Grid and axes
        if (this.options.showGrid) {
            const gridHelper = new this.THREE.GridHelper(20, 20);
            this.scene.add(gridHelper);
        }
        if (this.options.showAxes) {
            const axesHelper = new this.THREE.AxesHelper(5);
            this.scene.add(axesHelper);
        }

        // Controls
        this.controls = new this.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Resize handler
        window.addEventListener('resize', () => this.onResize());

        // Start animation
        this.animate();

        console.log('[Three3DRenderer] Initialized successfully');
    }

    animate() {
        if (!this.renderer || !this.scene || !this.camera) return;

        requestAnimationFrame(() => this.animate());

        if (this.controls) {
            this.controls.update();
        }

        this.renderer.render(this.scene, this.camera);
    }

    onResize() {
        if (!this.container || !this.camera || !this.renderer) return;

        const width = this.container.clientWidth;
        const height = this.container.clientHeight;

        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(width, height);
    }

    parseWKT(wkt) {
        const original = wkt;
        const wktUpper = wkt.trim().toUpperCase();

        // POINT (with or without Z)
        if (wktUpper.startsWith('POINT')) {
            const match = wkt.match(/POINT(?:\s+Z)?\s*\(\s*([^)]+)\s*\)/i);
            if (match) {
                const parts = match[1].trim().split(/\s+/).map(Number);
                return {
                    type: 'POINT',
                    coords: { x: parts[0], y: parts[1], z: parts[2] || 0 }
                };
            }
        }

        // LINESTRING (with or without Z)
        if (wktUpper.startsWith('LINESTRING')) {
            const match = wkt.match(/LINESTRING(?:\s+Z)?\s*\(([^)]+)\)/i);
            if (match) {
                const coords = match[1].split(',').map(point => {
                    const parts = point.trim().split(/\s+/).map(Number);
                    return { x: parts[0], y: parts[1], z: parts[2] || 0 };
                });
                return { type: 'LINESTRING', coords };
            }
        }

        // POLYGON (with or without Z) - first ring only (outer boundary)
        if (wktUpper.startsWith('POLYGON')) {
            const match = wkt.match(/POLYGON(?:\s+Z)?\s*\(\s*\(([^)]+)\)/i);
            if (match) {
                const coords = match[1].split(',').map(point => {
                    const parts = point.trim().split(/\s+/).map(Number);
                    return { x: parts[0], y: parts[1], z: parts[2] || 0 };
                });
                return { type: 'POLYGON', coords };
            }
        }

        // TRIANGLE (with or without Z)
        if (wktUpper.startsWith('TRIANGLE')) {
            const match = wkt.match(/TRIANGLE(?:\s+Z)?\s*\(\s*\(([^)]+)\)\s*\)/i);
            if (match) {
                const coords = match[1].split(',').map(point => {
                    const parts = point.trim().split(/\s+/).map(Number);
                    return { x: parts[0], y: parts[1], z: parts[2] || 0 };
                });
                return { type: 'TRIANGLE', coords };
            }
        }

        // POLYHEDRALSURFACE - extract all polygons
        if (wktUpper.startsWith('POLYHEDRALSURFACE')) {
            return { type: 'POLYHEDRALSURFACE', wkt: original };
        }

        // SOLID
        if (wktUpper.startsWith('SOLID')) {
            return { type: 'SOLID', wkt: original };
        }

        // TIN
        if (wktUpper.startsWith('TIN')) {
            return { type: 'TIN', wkt: original };
        }

        // GEOMETRYCOLLECTION - try to extract inner geometry
        if (wktUpper.startsWith('GEOMETRYCOLLECTION')) {
            // Try to extract TIN from GEOMETRYCOLLECTION Z (TIN Z (...))
            if (wktUpper.includes('TIN')) {
                const tinMatch = wkt.match(/TIN\s+Z?\s*\((.+)\)/is);
                if (tinMatch) {
                    return { type: 'TIN', wkt: 'TIN ' + tinMatch[0] };
                }
            }
            return { type: 'GEOMETRYCOLLECTION', wkt: original };
        }

        return { type: 'UNSUPPORTED', wkt: original };
    }

    extractPolygonsFromWKT(wkt) {
        const polygons = [];
        const polygonPattern = /\(\s*\(\s*([\d\s.,e+-]+)\s*\)\s*\)/g;
        let match;

        while ((match = polygonPattern.exec(wkt)) !== null) {
            const coords = [];
            const points = match[1].split(',');

            for (const point of points) {
                const parts = point.trim().split(/\s+/).filter(p => p.length > 0);
                if (parts.length >= 3) {
                    const x = parseFloat(parts[0]);
                    const y = parseFloat(parts[1]);
                    const z = parseFloat(parts[2]);
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        coords.push({ x, y, z });
                    }
                }
            }

            if (coords.length >= 3) {
                polygons.push(coords);
            }
        }

        return polygons;
    }

    createGeometryFromWKT(geom) {
        // POINT - small sphere
        if (geom.type === 'POINT') {
            const geometry = new this.THREE.SphereGeometry(0.15, 16, 16);
            geometry.translate(geom.coords.x, geom.coords.y, geom.coords.z);
            return geometry;
        }

        // LINESTRING - line
        if (geom.type === 'LINESTRING') {
            if (geom.coords.length < 2) {
                console.warn('[Three3DRenderer] LINESTRING has less than 2 points');
                return null;
            }
            const points = geom.coords.map(c => new this.THREE.Vector3(c.x, c.y, c.z));
            const geometry = new this.THREE.BufferGeometry().setFromPoints(points);
            return geometry;
        }

        // POLYGON or TRIANGLE - triangulated surface
        if (geom.type === 'POLYGON' || geom.type === 'TRIANGLE') {
            if (geom.coords.length < 3) {
                console.warn('[Three3DRenderer] Polygon has less than 3 coordinates');
                return null;
            }

            const vertices = [];
            const indices = [];

            geom.coords.forEach(coord => {
                vertices.push(coord.x, coord.y, coord.z);
            });

            // Fan triangulation from first vertex
            for (let i = 1; i < geom.coords.length - 1; i++) {
                indices.push(0, i, i + 1);
            }

            const geometry = new this.THREE.BufferGeometry();
            geometry.setAttribute('position', new this.THREE.Float32BufferAttribute(vertices, 3));
            geometry.setIndex(indices);
            geometry.computeVertexNormals();

            return geometry;
        }

        // POLYHEDRALSURFACE, SOLID, TIN - extract all polygons
        if (geom.type === 'POLYHEDRALSURFACE' || geom.type === 'SOLID' || geom.type === 'TIN') {
            const polygons = this.extractPolygonsFromWKT(geom.wkt);

            if (polygons.length === 0) {
                console.warn('[Three3DRenderer] No polygons extracted from', geom.type);
                return null;
            }

            const vertices = [];
            const indices = [];
            let vertexIndex = 0;

            for (const polygon of polygons) {
                if (polygon.length < 3) continue;

                const startIndex = vertexIndex;

                for (const coord of polygon) {
                    vertices.push(coord.x, coord.y, coord.z);
                    vertexIndex++;
                }

                // Fan triangulation from first vertex of polygon
                for (let i = 1; i < polygon.length - 1; i++) {
                    indices.push(startIndex, startIndex + i, startIndex + i + 1);
                }
            }

            if (vertices.length === 0) {
                console.warn('[Three3DRenderer] No vertices in', geom.type);
                return null;
            }

            const geometry = new this.THREE.BufferGeometry();
            geometry.setAttribute('position', new this.THREE.Float32BufferAttribute(vertices, 3));
            geometry.setIndex(indices);
            geometry.computeVertexNormals();

            return geometry;
        }

        console.warn('[Three3DRenderer] Unsupported geometry type:', geom.type);
        return null;
    }

    async addGeometry(wkt, opts = {}) {
        if (!this.initialized) {
            await this.initPromise;
        }

        if (!this.scene || !this.THREE) {
            console.warn('[Three3DRenderer] Not ready');
            return;
        }

        try {
            const geom = this.parseWKT(wkt);
            const geometry = this.createGeometryFromWKT(geom);

            if (!geometry) {
                console.warn('[Three3DRenderer] Could not create geometry from WKT');
                return;
            }

            const color = opts.color || this.colors[this.geometries.length % this.colors.length];
            let mesh;

            // LINESTRING uses Line instead of Mesh
            if (geom.type === 'LINESTRING') {
                const material = new this.THREE.LineBasicMaterial({
                    color: typeof color === 'string' ? parseInt(color.replace('#', '0x')) : color,
                    linewidth: 2
                });
                mesh = new this.THREE.Line(geometry, material);
            } else {
                // Regular mesh for surfaces
                const material = new this.THREE.MeshPhongMaterial({
                    color: typeof color === 'string' ? parseInt(color.replace('#', '0x')) : color,
                    side: this.THREE.DoubleSide,
                    transparent: true,
                    opacity: 0.8
                });

                mesh = new this.THREE.Mesh(geometry, material);

                // Add wireframe edges (not for lines or points)
                if (geom.type !== 'POINT') {
                    const edges = new this.THREE.EdgesGeometry(geometry);
                    const line = new this.THREE.LineSegments(
                        edges,
                        new this.THREE.LineBasicMaterial({ color: 0x000000 })
                    );
                    mesh.add(line);
                }
            }

            this.scene.add(mesh);
            this.geometries.push({ mesh, wkt });

            // Adjust camera to fit all geometries
            this.fitCameraToGeometries();

            console.log('[Three3DRenderer] Geometry added');
        } catch (error) {
            console.error('[Three3DRenderer] Error adding geometry:', error);
        }
    }

    fitCameraToGeometries() {
        if (this.geometries.length === 0) return;

        // Compute bounding box for all geometries
        const box = new this.THREE.Box3();
        this.geometries.forEach(({ mesh }) => {
            const meshBox = new this.THREE.Box3().setFromObject(mesh);
            box.union(meshBox);
        });

        if (box.isEmpty()) return;

        const center = new this.THREE.Vector3();
        box.getCenter(center);

        const size = new this.THREE.Vector3();
        box.getSize(size);
        const maxDim = Math.max(size.x, size.y, size.z);

        // Ensure minimum dimension for small objects (like single points)
        const distance = Math.max(maxDim * 2, 10);

        this.camera.position.set(distance, distance, distance);
        this.camera.lookAt(center);

        if (this.controls) {
            this.controls.target.copy(center);
            this.controls.update();
        }
    }

    setWireframe(enabled) {
        this.geometries.forEach(({ mesh }) => {
            if (mesh.material && mesh.material.wireframe !== undefined) {
                mesh.material.wireframe = enabled;
            }
        });
    }

    clear() {
        if (!this.scene) return;

        this.geometries.forEach(({ mesh }) => {
            this.scene.remove(mesh);
            if (mesh.geometry) mesh.geometry.dispose();
            if (mesh.material) mesh.material.dispose();
        });

        this.geometries = [];
    }
}
