using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class Scene
    {
        private UniformGrid _ug;

        public Scene(AABB bounds, int xdiv, int ydiv, int zdiv)
        {
            _ug = new UniformGrid(bounds, xdiv, ydiv, zdiv);
        }
        public void RegistTriangle(Triangle triangle)
        {
            _ug.RegistTriangle(triangle);
        }
        public IntersectionStatus RayTrace(Ray ray)
        {
            Triangle tri;
            float u, v, t;
            return this.RayTrace(ray, out tri, out u, out v, out t);
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3f p, out Vector3f n)
        {
            Triangle tri;
            float u, v, t;
            p = n = new Vector3f();

            IntersectionStatus stat = RayTrace(ray, out tri, out u, out v, out t);
            if (stat == IntersectionStatus.Intersect)
            {
                p = ray.Origin + t * ray.Direction;
                n = Vector3f.Normalize(tri.VA0.Normal + u * (tri.VA1.Normal - tri.VA0.Normal) + v * (tri.VA2.Normal - tri.VA0.Normal));
            }
            return stat;
        }
        public IntersectionStatus RayTrace(Ray ray, out Triangle triangle, out float u, out float v, out float t)
        {
            return _ug.RayTrace(ray, out triangle, out u, out v, out t);
        }

    };
    public class UniformGrid
    {
        public enum CubeFacetIdx
        {
            PositiveX,
            NegativeX,
            PositiveY,
            NegativeY,
            PositiveZ,
            NegativeZ,
        };

        public static CubeFacetIdx ToCubeFacetIdx(AABB.AABBFacet facet)
        {
            switch (facet)
            {
                case AABB.AABBFacet.PositiveX: return CubeFacetIdx.PositiveX;
                case AABB.AABBFacet.PositiveY: return CubeFacetIdx.PositiveY;
                case AABB.AABBFacet.PositiveZ: return CubeFacetIdx.PositiveZ;
                case AABB.AABBFacet.NegativeX: return CubeFacetIdx.NegativeX;
                case AABB.AABBFacet.NegativeY: return CubeFacetIdx.NegativeY;
                case AABB.AABBFacet.NegativeZ: return CubeFacetIdx.NegativeZ;
                default: throw new OverflowException();
            }
        }
     

        private UniformGridCell[, ,] _cell;
        private AABB _bounds;

        public int NumCells { get { return _cell.Length; } }
        public UniformGridCell[, ,] Cells { get { return _cell; } }
        public int NumCellsX { get { return _cell.GetLength(0); } }
        public int NumCellsY { get { return _cell.GetLength(1); } }
        public int NumCellsZ { get { return _cell.GetLength(2); } }

        public UniformGrid(AABB bounds, int xdiv, int ydiv, int zdiv)
        {
            _bounds = bounds;
            _cell = new UniformGridCell[xdiv, ydiv, zdiv];
            Vector3f cellHalfVector = new Vector3f(
                bounds.HalfVector.X * (float)1.0 / (float)xdiv,
                bounds.HalfVector.Y * (float)1.0 / (float)ydiv,
                bounds.HalfVector.Z * (float)1.0 / (float)zdiv);

            for (int i = 0; i < xdiv; i++)
            {
                for (int j = 0; j < ydiv; j++)
                {
                    for (int k = 0; k < zdiv; k++)
                    {
                        Vector3f center = new Vector3f(
                            (2.0f * (float)i + 1.0f) * cellHalfVector.X,
                            (2.0f * (float)j + 1.0f) * cellHalfVector.Y,
                            (2.0f * (float)k + 1.0f) * cellHalfVector.Z) + bounds.Min;
                        AABB aabb = new AABB(center, cellHalfVector);

                        _cell[i, j, k] = new UniformGridCell(aabb, i, j, k);
                    }
                }
            }
        }
        
        public void RegistTriangle(Triangle tri)
        {
            // 三角形が作るAABB
            AABB aabb = AABB.FromTriangel(tri);

            UniformGridCell[] cells = AabbOverlapCells(aabb);
            foreach (UniformGridCell cell in cells)
            {
                cell.RegistTriangle(tri);
            }
        }

        public IntersectionStatus RayTrace(Ray ray, out  Triangle triangle, out float u, out float v, out float t)
        {
            triangle = null;
            u = v = t = 0.0f;

            // Rayと衝突するすべてのCellを引っ張り出す。
            UniformGridCell[] cells = Traverse(ray);


            if (cells == null) return IntersectionStatus.Reject;

            for (int i = 0; i < cells.Length; i++)
            {
                if (cells[i].RayTrace(ray, out triangle, out u, out v, out t) == IntersectionStatus.Intersect)
                {
                    // 衝突点がこのセルの中にあるか判定
                    Vector3f p = ray.Origin + t * ray.Direction;
                    if (Intersections.PointAABBOverlap(p, cells[i].Bounds) == OverlapStatus.Disjoint) continue;
                    return IntersectionStatus.Intersect;

                }
            }

            return IntersectionStatus.Reject;

        }

        public UniformGridCell[] Traverse(Ray ray)
        {
            List<UniformGridCell> list = new List<UniformGridCell>();
            //{
            //    // debug
            //    foreach (UniformGridCell cell in _cell)
            //    {
            //        float t_min,t_max;
            //        AABB.AABBFacet nearFacet, farFacet;
            //        if (Intersections.RayAABBIntersect(ray.Origin, ray.Direction, cell.Bounds, out t_min, out t_max, out nearFacet, out farFacet) == IntersectionStatus.Intersect) list.Add(cell);
            //    }
            //    return list.ToArray();
            //    //
            //}

            {
                float t_min, t_max;
                AABB.AABBFacet nearFacet, farFacet;
                UniformGridCell cell;

                // Rayの始点のシーン内外判定
                if (Intersections.PointAABBOverlap(ray.Origin, _bounds) == OverlapStatus.Disjoint)
                {
                    // 外側のとき

                    // まずシーン全体のAABBのどこに当たるかを判定
                    if (Intersections.RayAABBIntersect(ray.Origin, ray.Direction, _bounds, out t_min, out  t_max, out nearFacet, out farFacet) == IntersectionStatus.Reject) return null;
                    else
                    {
                        // 衝突点（近い方）
                        Vector3f pos = ray.Origin + t_min * ray.Direction;

                        // 衝突点を含むCellを抽出
                        cell = PointOnBoundsOverlapCell(pos, ToCubeFacetIdx(nearFacet));
                    }
                }
                else
                {
                    // 内側のとき

                    // 始点を含むCellを抽出
                    cell = PointOverlapCell(ray.Origin);
                }

                // 外側の衝突面
                CubeFacetIdx idx;

                while (cell != null)
                {
                    // 移動したセルと衝突判定
                    if (Intersections.RayAABBIntersect(ray.Origin, ray.Direction, cell.Bounds, out t_min, out t_max, out nearFacet, out farFacet) == IntersectionStatus.Reject)
                    {
                        // なぜここを通過するやつがいる？
                        break;
                    }
                    list.Add(cell);

                    idx = ToCubeFacetIdx(farFacet);

                    //隣のセルに移動
                    cell = NeighborCell(cell, idx);
                }

                return list.ToArray();
            }
        }

        private UniformGridCell PointOnBoundsOverlapCell(Vector3f p, CubeFacetIdx idx)
        {
            // シーンのUGの大きさ
            float unit_x = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.X;
            float unit_y = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Y;
            float unit_z = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Z;

            // シーンの最小
            Vector3f min = _bounds.Min;

            // 最小からの距離
            Vector3f len = p - min;

            int id_x = (int)Math.Floor(len.X / unit_x);
            int id_y = (int)Math.Floor(len.Y / unit_y);
            int id_z = (int)Math.Floor(len.Z / unit_z);

            if (idx == CubeFacetIdx.PositiveX) id_x = NumCellsX - 1;
            else if (idx == CubeFacetIdx.NegativeX) id_x = 0;
            else if (idx == CubeFacetIdx.PositiveY) id_y = NumCellsY - 1;
            else if (idx == CubeFacetIdx.NegativeY) id_y = 0;
            else if (idx == CubeFacetIdx.PositiveZ) id_z = NumCellsZ - 1;
            else if (idx == CubeFacetIdx.NegativeZ) id_z = 0;

            return _cell[id_x, id_y, id_z];
        }

        /// <summary>
        /// 隣のセルを返す
        /// </summary>
        /// <param name="cell"></param>
        /// <param name="idx"></param>
        /// <returns></returns>
        private UniformGridCell NeighborCell(UniformGridCell cell, CubeFacetIdx idx)
        {
            int x_idx = cell.X;
            int y_idx = cell.Y;
            int z_idx = cell.Z;

            if (idx == CubeFacetIdx.PositiveX) x_idx++;
            else if (idx == CubeFacetIdx.NegativeX) x_idx--;
            else if (idx == CubeFacetIdx.PositiveY) y_idx++;
            else if (idx == CubeFacetIdx.NegativeY) y_idx--;
            else if (idx == CubeFacetIdx.PositiveZ) z_idx++;
            else if (idx == CubeFacetIdx.NegativeZ) z_idx--;

            if (x_idx >= NumCellsX || x_idx < 0) return null;
            else if (y_idx >= NumCellsY || y_idx < 0) return null;
            else if (z_idx >= NumCellsZ || z_idx < 0) return null;

            return _cell[x_idx, y_idx, z_idx];

        }

        private UniformGridCell PointOverlapCell(Vector3f p)
        {
            // シーンのUGの大きさ
            float unit_x = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.X;
            float unit_y = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Y;
            float unit_z = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Z;

            // シーンの最小
            Vector3f min = _bounds.Min;

            // 最小からの距離
            Vector3f len = p - min;

            int id_x = (int)Math.Floor(len.X / unit_x);
            int id_y = (int)Math.Floor(len.Y / unit_y);
            int id_z = (int)Math.Floor(len.Z / unit_z);

            return _cell[id_x, id_y, id_z];
        }

        /// <summary>
        /// AABBとオーバーラップするCELLを返却する
        /// </summary>
        /// <param name="aabb"></param>
        /// <returns></returns>
        private UniformGridCell[] AabbOverlapCells(AABB aabb)
        {
            List<UniformGridCell> list = new List<UniformGridCell>();

            // 指定のAABBの最小
            float min_x = aabb.Min.X;
            float min_y = aabb.Min.Y;
            float min_z = aabb.Min.Z;

            // シーンのUGの大きさ
            float unit_x = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.X;
            float unit_y = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Y;
            float unit_z = 2.0f * _cell[0, 0, 0].Bounds.HalfVector.Z;

            // シーンの最小
            float bounds_x = _bounds.Min.X;
            float bounds_y = _bounds.Min.Y;
            float bounds_z = _bounds.Min.Z;

            // 最小値から指定のAABBの最小までの距離
            float len_x = min_x - bounds_x;
            float len_y = min_y - bounds_y;
            float len_z = min_z - bounds_z;

            // そこにCELLはいくつ入る
            int num_cell_min_x = (int)Math.Floor(len_x / unit_x);
            int num_cell_max_x = (int)Math.Ceiling((aabb.Max.X - _bounds.Min.X) / unit_x);
            int num_cell_min_y = (int)Math.Floor(len_y / unit_y);
            int num_cell_max_y = (int)Math.Ceiling((aabb.Max.Y - _bounds.Min.Y) / unit_y);
            int num_cell_min_z = (int)Math.Floor(len_z / unit_z);
            int num_cell_max_z = (int)Math.Ceiling((aabb.Max.Z - _bounds.Min.Z) / unit_z);

            int num_div_x = this.NumCellsX;
            int num_div_y = this.NumCellsY;
            int num_div_z = this.NumCellsZ;

            for (int i = num_cell_min_x; i < num_cell_max_x; i++)
            {
                for (int j = num_cell_min_y; j < num_cell_max_y; j++)
                {
                    for (int k = num_cell_min_z; k < num_cell_max_z; k++)
                    {
                        list.Add(_cell[i, j, k]);
                    }
                }
            }

            return list.ToArray();
        }

    };

    public class UniformGridCell
    {
        private AABB _bounds;
        private List<Triangle> _triangles = new List<Triangle>();
        private int _xid, _yid, _zid;

        public AABB Bounds { get { return _bounds; } }
        public int X { get { return _xid; } }
        public int Y { get { return _yid; } }
        public int Z { get { return _zid; } }

        public UniformGridCell(AABB bounds, int xid, int yid, int zid)
        {
            _bounds = bounds;
            _xid = xid;
            _yid = yid;
            _zid = zid;
        }

        public void RegistTriangle(Triangle tri)
        {
            //if (Intersections.TriangleAABBOverlap(this.Bounds, tri) == OverlapStatus.Overlap) _triangles.Add(tri);
            _triangles.Add(tri);
        }

        internal IntersectionStatus RayTrace(Ray ray, out Triangle triangle, out float u, out float v, out float t)
        {
            triangle = null;
            u = v = t = 0.0f;

            float t_min = float.MaxValue;
            int id = -1;

            for (int i = 0; i < _triangles.Count; i++)
            {
                float temp_u, temp_v, temp_t;
                if (Intersections.RayTriangleIntersect(ray, _triangles[i], out temp_u, out temp_v, out temp_t) == IntersectionStatus.Intersect)
                {
                    // 始点より後ろでぶつかってたら無視
                    if (temp_t < 0.0) continue;

                    if (temp_t < t_min)
                    {
                        t_min = temp_t;
                        id = i;

                        t = temp_t;
                        u = temp_u;
                        v = temp_v;
                    }
                }
            }

            if (id != -1)
            {
                triangle = _triangles[id];
                return IntersectionStatus.Intersect;
            }
            else return IntersectionStatus.Reject;
        }
    }
}
