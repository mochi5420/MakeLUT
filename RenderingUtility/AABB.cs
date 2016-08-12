using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class AABB
    {
        public enum AABBFacet
        {
            PositiveX,
            PositiveY,
            PositiveZ,
            NegativeX,
            NegativeY,
            NegativeZ
        };


        private Vector3f _center;
        private Vector3f _half;
        private static int[] _vertexIndices = new int[] { 2, 0, 1, 1, 3, 2, 4, 0, 2, 2, 6, 4,
            5, 1, 0, 0, 4, 5, 7, 3, 1, 1, 5, 7, 6, 2, 3, 3, 7, 6, 4, 6, 7, 7, 5, 4 };
        private static int[] _vertexIndicesLines = new int[]{
            0,1,1,3,3,2,2,0,4,5,5,7,7,6,6,4,0,4,1,5,2,6,3,7
        };

        public static int[] VertexIndices
        {
            get { return AABB._vertexIndices; }
        }
        public static int[] LineListVertexIndices
        {
            get { return AABB._vertexIndicesLines; }
        }
        public AABB(float xMax, float xMin, float yMax, float yMin, float zMax, float zMin)
        {
            _center = 0.5f * new Vector3f(xMax + xMin, yMax + yMin, zMax + zMin);
            _half = new Vector3f(xMax, yMax, zMax) - _center;
        }
        public AABB(Vector3f center, Vector3f halfvector)
        {
            _center = center;
            _half = halfvector;
        }
        public Vector3f Center
        {
            set { _center = value; }
            get { return _center; }
        }
        public Vector3f HalfVector
        {
            set { _half = value; }
            get { return _half; }
        }
        public Vector3f Vertex(int idx)
        {
            switch (idx)
            {
                case 0:
                    return _center + new Vector3f(-_half.X, +_half.Y, +_half.Z);
                case 1:
                    return _center + new Vector3f(+_half.X, +_half.Y, +_half.Z);
                case 2:
                    return _center + new Vector3f(-_half.X, +_half.Y, -_half.Z);
                case 3:
                    return _center + new Vector3f(+_half.X, +_half.Y, -_half.Z);
                case 4:
                    return _center + new Vector3f(-_half.X, -_half.Y, +_half.Z);
                case 5:
                    return _center + new Vector3f(+_half.X, -_half.Y, +_half.Z);
                case 6:
                    return _center + new Vector3f(-_half.X, -_half.Y, -_half.Z);
                case 7:
                    return _center + new Vector3f(+_half.X, -_half.Y, -_half.Z);
                default:
                    throw new Exception("Unsupported Index.");
            }
        }
        public Vector3f[] TriangleStripVertices
        {
            get
            {
                Vector3f[] p = new Vector3f[36];

                // Positive-X
                p[0] = Vertex(5); p[1] = Vertex(7); p[2] = Vertex(3);
                p[3] = Vertex(5); p[4] = Vertex(3); p[5] = Vertex(1);

                // Negative-X
                p[6] = Vertex(0); p[7] = Vertex(2); p[8] = Vertex(6);
                p[9] = Vertex(0); p[10] = Vertex(6); p[11] = Vertex(4);

                // Positive-Y
                p[12] = Vertex(0); p[13] = Vertex(1); p[14] = Vertex(2);
                p[15] = Vertex(1); p[16] = Vertex(3); p[17] = Vertex(2);

                // Negative-Y
                p[18] = Vertex(6); p[19] = Vertex(5); p[20] = Vertex(4);
                p[21] = Vertex(6); p[22] = Vertex(7); p[23] = Vertex(5);

                // Positive-Z
                p[24] = Vertex(0); p[25] = Vertex(5); p[26] = Vertex(1);
                p[27] = Vertex(0); p[28] = Vertex(4); p[29] = Vertex(5);

                // Negative-Z
                p[30] = Vertex(7); p[31] = Vertex(2); p[32] = Vertex(3);
                p[33] = Vertex(7); p[34] = Vertex(6); p[35] = Vertex(2);

                return p;
            }
        }

        public float Right
        {
            get { return _center.X + _half.X; }
        }
        public float Left
        {
            get { return _center.X - _half.X; }
        }
        public float Front
        {
            get { return _center.Y + _half.Y; }
        }
        public float Back
        {
            get { return _center.Y - _half.Y; }
        }

        public float Top
        {
            get { return _center.Z + _half.Z; }
        }
        public float Bottom
        {
            get { return _center.Z - _half.Z; }
        }

        public Vector3f Min
        {
            get { return _center - _half; }
        }
        public Vector3f Max
        {
            get { return _center + _half; }
        }

        static public AABB FromTriangel(Triangle tri)
        {
            Vector3f max = new Vector3f(
                Math.Max(tri.V0.X, Math.Max(tri.V1.X, tri.V2.X)),
                Math.Max(tri.V0.Y, Math.Max(tri.V1.Y, tri.V2.Y)),
                Math.Max(tri.V0.Z, Math.Max(tri.V1.Z, tri.V2.Z)));
            Vector3f min = new Vector3f(
                Math.Min(tri.V0.X, Math.Min(tri.V1.X, tri.V2.X)),
                Math.Min(tri.V0.Y, Math.Min(tri.V1.Y, tri.V2.Y)),
                Math.Min(tri.V0.Z, Math.Min(tri.V1.Z, tri.V2.Z)));

            Vector3f center = 0.5f * (max + min);


            return new AABB(center, max - center);
        }


    }
}
