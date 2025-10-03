# SFCGAL WebAssembly - Known Limitations

## Memory Constraints

WebAssembly has memory limitations that affect complex 3D operations:

- **Maximum memory**: ~2GB (4GB on 64-bit browsers with specific flags)
- **Initial memory**: 256MB (configurable)
- **Memory growth**: Enabled, but can be slow for large allocations

## 3D Boolean Operations (intersection3D, union3D, difference3D)

### ⚠️ IMPORTANT: Use 2D algorithms for planar geometries

CGAL's 3D boolean operations use **exact arithmetic with Nef polyhedra**, which can fail with:
- **"Undecidable conversion of CGAL::Uncertain<T>"**: When geometric predicates cannot be determined with certainty
- This commonly occurs with **planar geometries** (all points have same Z coordinate)

**Solution**: For planar geometries, use the 2D algorithms instead:
- ✅ Use `intersection()` instead of `intersection3D()`
- ✅ Use `union()` instead of `union3D()`
- ✅ Use `difference()` instead of `difference3D()`

The 2D algorithms work perfectly with `POLYGON Z` geometries as long as they're planar.

### Working Well ✅
- **POLYGON Z**: Simple 3D polygons work efficiently
  ```
  POLYGON Z ((0 0 0, 4 0 0, 4 3 0, 0 3 0, 0 0 0))
  ```
- **Small POLYHEDRALSURFACE**: Simple polyhedral surfaces with few faces
- **Extruded simple shapes**: Basic extrusions with low complexity

### Memory-Intensive ⚠️
- **SOLID geometries**: Complex solids require significant memory
  - Operations can require 1-2GB+ for moderately complex solids
  - May cause "Cannot allocate memory" errors
  - May cause "index out of bounds" errors with invalid topology

### Recommended Approach

For WebAssembly demonstrations:
1. Use **POLYGON Z** for 3D boolean operations examples
2. Keep geometries simple (< 100 vertices per geometry)
3. Use **POLYHEDRALSURFACE** only for simple shapes (< 10 faces)
4. Avoid **SOLID** operations in browser examples

For production use with complex 3D operations:
- Use SFCGAL in **Node.js** with `--max-old-space-size` flag
- Use SFCGAL **native binaries** (C++ or PostGIS)
- Consider **server-side processing** for complex operations

## Error Messages

### "Cannot allocate memory (size=...)"
- The operation requires more memory than available
- Solution: Simplify geometries or use simpler operation

### "index out of bounds"
- SFCGAL detected invalid geometry topology
- Solution: Validate geometry with `isValid()` before operation
- May indicate WKT parsing error

### "Aborted()"
- Critical error in CGAL/SFCGAL algorithm
- Usually due to memory exhaustion or invalid input
- Check browser console for detailed error messages

## Best Practices

1. **Always validate geometries** before complex operations:
   ```javascript
   if (!sfcgal.isValid(wkt)) {
       console.error("Invalid geometry");
   }
   ```

2. **Use appropriate precision** in WKT (1 decimal place is often sufficient):
   ```javascript
   result.asText(1)  // Good
   result.asText(15) // May increase memory usage
   ```

3. **Catch errors gracefully**:
   ```javascript
   try {
       const result = sfcgal.intersection3D(wkt1, wkt2);
   } catch (error) {
       console.error("Operation failed:", error.message);
   }
   ```

4. **Monitor browser console** for warnings from SFCGAL binding

## Performance Tips

- 2D operations are generally faster and use less memory
- `area()`, `length()`, `distance()` are lightweight
- Boolean operations (union, intersection, difference) are more expensive
- 3D operations are significantly more expensive than 2D equivalents
- `convexHull()` is relatively fast even for complex geometries
