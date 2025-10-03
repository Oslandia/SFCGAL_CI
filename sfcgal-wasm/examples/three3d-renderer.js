/**
 * Three.js 3D Renderer for SFCGAL Geometries
 * Renders WKT 3D geometries with Three.js (solid/wireframe modes)
 * Using ES6 modules (Three.js r160+)
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

export class Three3DRenderer {
    constructor(containerId, options = {}) {
        this.container = document.getElementById(containerId);
        if (!this.container) {
            throw new Error(`Container element #${containerId} not found`);
        }

        this.options = {
            wireframe: options.wireframe || false,
            solidColor: options.solidColor || 0x667eea,
            wireframeColor: options.wireframeColor || 0x333333,
            gridSize: options.gridSize || 20,
            showGrid: options.showGrid !== false,
            showAxes: options.showAxes !== false,
            backgroundColor: options.backgroundColor || 0xf5f5f5,
            cameraDistance: options.cameraDistance || 10
        };

        this.initScene();
        this.geometries = [];
        this.animate();
    }

    initScene() {
        // Scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(this.options.backgroundColor);

        // Camera
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(60, width / height, 0.1, 1000);
        this.camera.position.set(5, 5, 5);

        // Renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(width, height);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        // Lights
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        this.scene.add(ambientLight);

        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(10, 10, 10);
        this.scene.add(directionalLight);

        // Grid
        if (this.options.showGrid) {
            const gridHelper = new THREE.GridHelper(this.options.gridSize, this.options.gridSize);
            gridHelper.material.color.setHex(0xcccccc);
            this.scene.add(gridHelper);
        }

        // Axes
        if (this.options.showAxes) {
            const axesHelper = new THREE.AxesHelper(this.options.gridSize / 2);
            this.scene.add(axesHelper);
        }

        // Controls
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Handle resize
        window.addEventListener('resize', () => this.onResize());
    }

    onResize() {
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;

        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(width, height);
    }

    clear() {
        // Remove all geometry meshes
        for (const mesh of this.geometries) {
            this.scene.remove(mesh);
            if (mesh.geometry) mesh.geometry.dispose();
            if (mesh.material) mesh.material.dispose();
        }
        this.geometries = [];
    }

    setWireframe(enabled) {
        this.options.wireframe = enabled;
        // Update existing geometries
        for (const mesh of this.geometries) {
            if (mesh.material) {
                mesh.material.wireframe = enabled;
            }
        }
    }

    addGeometry(wkt, options = {}) {
        const mesh = this.wktToMesh(wkt, options);
        if (mesh) {
            this.scene.add(mesh);
            this.geometries.push(mesh);
            this.fitCameraToGeometries();
        }
    }

    wktToMesh(wkt, options = {}) {
        const upperWKT = wkt.toUpperCase();
        const wireframe = options.wireframe !== undefined ? options.wireframe : this.options.wireframe;
        const color = options.color || this.options.solidColor;

        if (upperWKT.startsWith('GEOMETRYCOLLECTION')) {
            return this.createGeometryCollectionMesh(wkt, color, wireframe);
        } else if (upperWKT.includes('TIN')) {
            return this.createTinMesh(wkt, color, wireframe);
        } else if (upperWKT.includes('POINT')) {
            return this.createPointMesh(wkt, color);
        } else if (upperWKT.includes('LINESTRING')) {
            return this.createLineMesh(wkt, color);
        } else if (upperWKT.includes('TRIANGLE') && !upperWKT.includes('TRIANGULATEDSURFACE')) {
            return this.createTriangleMesh(wkt, color, wireframe);
        } else if (upperWKT.includes('POLYGON') && !upperWKT.includes('POLYHEDRALSURFACE')) {
            return this.createPolygonMesh(wkt, color, wireframe);
        } else if (upperWKT.includes('SOLID') || upperWKT.includes('POLYHEDRALSURFACE')) {
            return this.createSolidMesh(wkt, color, wireframe);
        }

        return null;
    }

    extractCoords(wkt) {
        const coords = [];
        // Match X Y Z coordinates (handles optional Z)
        const regex = /-?\d+\.?\d*\s+-?\d+\.?\d*(?:\s+-?\d+\.?\d*)?/g;
        const matches = wkt.match(regex);

        if (matches) {
            for (const match of matches) {
                const parts = match.trim().split(/\s+/).map(Number);
                coords.push({
                    x: parts[0] || 0,
                    y: parts[1] || 0,
                    z: parts[2] || 0
                });
            }
        }

        return coords;
    }

    createPointMesh(wkt, color) {
        const coords = this.extractCoords(wkt);
        if (coords.length === 0) return null;

        const geometry = new THREE.SphereGeometry(0.2, 16, 16);
        const material = new THREE.MeshPhongMaterial({ color });
        const mesh = new THREE.Mesh(geometry, material);

        mesh.position.set(coords[0].x, coords[0].z, -coords[0].y); // Y/Z swap for Three.js
        return mesh;
    }

    createLineMesh(wkt, color) {
        const coords = this.extractCoords(wkt);
        if (coords.length < 2) return null;

        const points = coords.map(c => new THREE.Vector3(c.x, c.z, -c.y));
        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({ color });

        return new THREE.Line(geometry, material);
    }

    createTriangleMesh(wkt, color, wireframe) {
        const coords = this.extractCoords(wkt);
        if (coords.length < 3) return null;

        // TRIANGLE can be 2D or 3D - check if it has Z coordinates
        const is3D = coords.some(c => c.z !== 0);

        const geometry = new THREE.BufferGeometry();
        const vertices = [];
        const indices = [];

        // Add vertices with Y/Z swap for Three.js coordinate system
        for (let i = 0; i < Math.min(3, coords.length); i++) {
            if (is3D) {
                vertices.push(coords[i].x, coords[i].z, -coords[i].y);
            } else {
                // 2D triangle - use x, y, 0 for z
                vertices.push(coords[i].x, 0, -coords[i].y);
            }
        }

        // Single triangle (3 vertices)
        indices.push(0, 1, 2);

        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color,
            wireframe,
            side: THREE.DoubleSide
        });

        return new THREE.Mesh(geometry, material);
    }

    createPolygonMesh(wkt, color, wireframe) {
        const coords = this.extractCoords(wkt);
        if (coords.length < 3) return null;

        // Check if it's a 3D polygon (has Z coordinates)
        const is3D = coords.some(c => c.z !== 0);

        if (is3D) {
            // Create a 3D polygon using BufferGeometry
            const geometry = new THREE.BufferGeometry();
            const vertices = [];
            const indices = [];

            // Add vertices with Y/Z swap for Three.js coordinate system
            for (const coord of coords) {
                vertices.push(coord.x, coord.z, -coord.y);
            }

            // Triangulate polygon (simple fan triangulation from first vertex)
            for (let i = 1; i < coords.length - 2; i++) {
                indices.push(0, i, i + 1);
            }

            geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
            geometry.setIndex(indices);
            geometry.computeVertexNormals();

            const material = new THREE.MeshPhongMaterial({
                color,
                wireframe,
                side: THREE.DoubleSide
            });

            return new THREE.Mesh(geometry, material);
        } else {
            // 2D polygon - use shape geometry
            const shape = new THREE.Shape();
            shape.moveTo(coords[0].x, coords[0].y);
            for (let i = 1; i < coords.length; i++) {
                shape.lineTo(coords[i].x, coords[i].y);
            }

            const geometry = new THREE.ShapeGeometry(shape);
            const material = new THREE.MeshPhongMaterial({
                color,
                wireframe,
                side: THREE.DoubleSide
            });

            const mesh = new THREE.Mesh(geometry, material);
            mesh.rotation.x = -Math.PI / 2; // Lay flat

            return mesh;
        }
    }

    createSolidMesh(wkt, color, wireframe) {
        // Parse POLYHEDRALSURFACE - extract individual polygons
        const polygonRegex = /\(\(([^)]+)\)\)/g;
        const polygons = [];
        let match;

        while ((match = polygonRegex.exec(wkt)) !== null) {
            const coordStr = match[1];
            const coords = [];
            const coordRegex = /-?\d+\.?\d*\s+-?\d+\.?\d*(?:\s+-?\d+\.?\d*)?/g;
            let coordMatch;

            while ((coordMatch = coordRegex.exec(coordStr)) !== null) {
                const parts = coordMatch[0].trim().split(/\s+/).map(Number);
                coords.push({
                    x: parts[0] || 0,
                    y: parts[1] || 0,
                    z: parts[2] || 0
                });
            }

            if (coords.length >= 3) {
                polygons.push(coords);
            }
        }

        if (polygons.length === 0) return null;

        // Build geometry from polygons
        const geometry = new THREE.BufferGeometry();
        const vertices = [];
        const indices = [];
        let vertexIndex = 0;

        for (const polygon of polygons) {
            // Triangulate polygon (simple fan triangulation from first vertex)
            const startIndex = vertexIndex;

            // Add all vertices of this polygon
            for (const coord of polygon) {
                vertices.push(coord.x, coord.z, -coord.y); // Y/Z swap for Three.js
                vertexIndex++;
            }

            // Create triangles (fan from first vertex)
            for (let i = 1; i < polygon.length - 1; i++) {
                indices.push(startIndex, startIndex + i, startIndex + i + 1);
            }
        }

        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color,
            wireframe,
            side: THREE.DoubleSide
        });

        return new THREE.Mesh(geometry, material);
    }

    createTinMesh(wkt, color, wireframe) {
        // TIN (Triangulated Irregular Network) - collection of triangles
        // Format: TIN Z (((x1 y1 z1,x2 y2 z2,x3 y3 z3,x1 y1 z1)),((...),...))

        // Remove TIN Z ( and trailing )
        let content = wkt.replace(/^TIN\s*Z?\s*\(/i, '');
        content = content.replace(/\)$/, '');

        // Parse triangles - each triangle is ((coords))
        const triangles = [];
        let depth = 0;
        let start = -1;

        for (let i = 0; i < content.length; i++) {
            const char = content[i];

            if (char === '(') {
                if (depth === 0) {
                    start = i + 1; // Start after the first '('
                }
                depth++;
            } else if (char === ')') {
                depth--;
                if (depth === 0 && start !== -1) {
                    // End of a triangle - extract coordinates between (( and ))
                    const triangleStr = content.substring(start, i);
                    // Remove the inner '(' if present
                    const coordStr = triangleStr.replace(/^\s*\(/, '').replace(/\)\s*$/, '');

                    const coords = [];
                    const coordRegex = /-?\d+\.?\d*\s+-?\d+\.?\d*(?:\s+-?\d+\.?\d*)?/g;
                    let coordMatch;

                    while ((coordMatch = coordRegex.exec(coordStr)) !== null) {
                        const parts = coordMatch[0].trim().split(/\s+/).map(Number);
                        coords.push({
                            x: parts[0] || 0,
                            y: parts[1] || 0,
                            z: parts[2] || 0
                        });
                    }

                    // Each triangle should have 3 or 4 points (last = first if 4)
                    if (coords.length >= 3) {
                        triangles.push(coords.slice(0, 3)); // Use first 3 points
                    }

                    start = -1;
                }
            }
        }

        if (triangles.length === 0) {
            console.warn('No triangles found in TIN:', wkt);
            return null;
        }

        console.log(`TIN: Found ${triangles.length} triangles`);

        // Build geometry from triangles
        const geometry = new THREE.BufferGeometry();
        const vertices = [];
        const indices = [];
        let vertexIndex = 0;

        for (const triangle of triangles) {
            const startIndex = vertexIndex;

            // Add vertices with Y/Z swap for Three.js
            for (const coord of triangle) {
                vertices.push(coord.x, coord.z, -coord.y);
                vertexIndex++;
            }

            // Add triangle indices
            indices.push(startIndex, startIndex + 1, startIndex + 2);
        }

        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();

        const material = new THREE.MeshPhongMaterial({
            color,
            wireframe,
            side: THREE.DoubleSide
        });

        return new THREE.Mesh(geometry, material);
    }

    createGeometryCollectionMesh(wkt, color, wireframe) {
        // GEOMETRYCOLLECTION can contain any mix of geometries
        // We need to properly parse nested parentheses
        const group = new THREE.Group();

        // Remove the GEOMETRYCOLLECTION wrapper
        let content = wkt.replace(/^GEOMETRYCOLLECTION\s*Z?\s*\(/i, '');
        content = content.replace(/\)$/, '');

        // Parse geometries by matching balanced parentheses
        const geometries = this.extractGeometriesFromCollection(content);

        console.log(`GEOMETRYCOLLECTION: Found ${geometries.length} sub-geometries`);

        for (const geomWkt of geometries) {
            console.log('Parsing sub-geometry:', geomWkt.substring(0, 50) + '...');
            const mesh = this.wktToMesh(geomWkt.trim(), { color, wireframe });
            if (mesh) {
                if (mesh instanceof THREE.Group) {
                    // If it's a group, add all children
                    for (const child of mesh.children) {
                        group.add(child);
                    }
                } else {
                    group.add(mesh);
                }
            } else {
                console.warn('Failed to create mesh for:', geomWkt.substring(0, 50) + '...');
            }
        }

        console.log(`GEOMETRYCOLLECTION: Added ${group.children.length} meshes to scene`);
        return group.children.length > 0 ? group : null;
    }

    extractGeometriesFromCollection(content) {
        // Extract individual geometries from a GEOMETRYCOLLECTION content
        const geometries = [];
        let depth = 0;
        let start = 0;
        let inGeometry = false;

        for (let i = 0; i < content.length; i++) {
            const char = content[i];

            if (char === '(') {
                if (depth === 0 && !inGeometry) {
                    // Find geometry type before this opening parenthesis
                    const before = content.substring(0, i).trim();
                    const typeMatch = before.match(/[A-Z]+\s*Z?\s*$/);
                    if (typeMatch) {
                        start = before.lastIndexOf(typeMatch[0]);
                        inGeometry = true;
                    }
                }
                depth++;
            } else if (char === ')') {
                depth--;
                if (depth === 0 && inGeometry) {
                    // End of a geometry
                    geometries.push(content.substring(start, i + 1));
                    inGeometry = false;
                }
            }
        }

        return geometries;
    }

    fitCameraToGeometries() {
        if (this.geometries.length === 0) return;

        const box = new THREE.Box3();
        for (const mesh of this.geometries) {
            box.expandByObject(mesh);
        }

        const center = box.getCenter(new THREE.Vector3());
        const size = box.getSize(new THREE.Vector3());
        const maxDim = Math.max(size.x, size.y, size.z);
        const fov = this.camera.fov * (Math.PI / 180);
        const cameraDistance = Math.abs(maxDim / Math.sin(fov / 2)) * 1.5;

        this.camera.position.set(
            center.x + cameraDistance,
            center.y + cameraDistance,
            center.z + cameraDistance
        );

        this.controls.target.copy(center);
        this.controls.update();
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }
}
