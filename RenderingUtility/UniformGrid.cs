using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class UniformGrid
    {
        private AABB _bounds;
        private UniformGridCell[, ,] _cells;

        public UniformGridCell[, ,] Cells
        {
            get { return _cells; }
            set { _cells = value; }
        }

        public AABB Bounds { get { return _bounds; } }


        public UniformGrid(AABB bounds, int xdiv, int ydiv, int zdiv)
        {
            _bounds = bounds;

            double unit_x = 1.0 / (double)xdiv * 2.0 * bounds.HalfVector.X;
            double unit_y = 1.0 / (double)ydiv * 2.0 * bounds.HalfVector.Y;
            double unit_z = 1.0 / (double)zdiv * 2.0 * bounds.HalfVector.Z;

            _cells = new UniformGridCell[xdiv, ydiv, zdiv];

            for (int i = 0; i < xdiv; i++)
            {
                for (int j = 0; j < ydiv; j++)
                {
                    for (int k = 0; k < zdiv; k++)
                    {
                        Vector3d half = 0.5 * new Vector3d(unit_x, unit_y, unit_z);
                        Vector3d center = bounds.Center - bounds.HalfVector
                            + new Vector3d((2.0 * (double)i + 0.5) * half.X, (2.0 * (double)j + 0.5) * half.Y, (2.0 * (double)k + 0.5) * half.Z);
                        _cells[i, j, k] = new UniformGridCell(new AABB(center, half));
                    }
                }
            }
        }


        private List<UniformGridCell> FindCells(Ray ray)
        {
            List<UniformGridCell> list = new List<UniformGridCell>();

            // まず全体とぶつける
            double t_near, t_far;
            AABB.AABBFacet facet_near, facet_far;
            if (Intersections.RayAABBIntersect(ray.Origin, ray.Direction, _bounds, out t_near, out t_far, out facet_near, out facet_far) == IntersectionStatus.Intersect)
            {
                // 衝突点
                Vector3d p = ray.Origin + t_near * ray.Direction;



                return list;
            }
            else
            {
                return null;
            }
        }
    }

    public class UniformGridCell
    {
        private AABB _bounds;

        public AABB Bounds
        {
            get { return _bounds; }
            set { _bounds = value; }
        }
        private List<Triangle> _triangles = new List<Triangle>();

        public UniformGridCell(AABB bounds)
        {
            _bounds = bounds;
        }

        private void RegistTriangle(Triangle triangle)
        {
            throw new NotImplementedException();
        }

    }
}
