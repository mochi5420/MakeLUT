using System;
using System.Collections.Generic;
using System.Text;
//using System.Collections;
using Utility;
using System.Xml.Serialization;

namespace RenderingUtility
{
    public enum CullMode
    {
        CullNone,
        Clockwise,
        CounterClockwise,
    };

    public class OctreeScene
    {

        private Octree _octree;


        // ‚RŠpŒ`
        private List<Triangle> _triangles = new List<Triangle>();

        public void SetupOctree(AABB aabb, int level)
        {
            _octree = new Octree(aabb);

            // ‹óŠÔ‚Í‚ ‚ç‚©‚¶‚ß•ªŠ„‚µ‚Ä‚¨‚­B
            for (int i = 0; i < level; i++) _octree.SplitAllLeaves();
        }
        public void SetupOctree(int level)
        {
            // ‰¼’u‚«‚³‚ê‚Ä‚¢‚éTriangle‚©‚ç©“®‚ÅAABB‚ğŒvZ
            double min_x = double.MaxValue;
            double min_y = double.MaxValue;
            double min_z = double.MaxValue;
 
            double max_x = double.MinValue;
            double max_y = double.MinValue;
            double max_z = double.MinValue;

            for (int i = 0; i < _triangles.Count; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double x = _triangles[i].V[j].X;
                    double y = _triangles[i].V[j].Y;
                    double z = _triangles[i].V[j].Z;

                    if (min_x > x) min_x = x;
                    if (min_y > y) min_y = y;
                    if (min_z > z) min_z = z;

                    if (max_x < x) max_x = x;
                    if (max_y < y) max_y = y;
                    if (max_z < z) max_z = z;
                }
            }

            double epsilon = 0.001; //‹C‚¿‘å‚«‚ß‚Éì‚Á‚Ä‚¨‚­B
            AABB aabb = new AABB(max_x + epsilon, min_x - epsilon, max_y + epsilon, min_y - epsilon, max_z + epsilon, min_z - epsilon);
            _octree = new Octree(aabb);
            for (int i = 0; i < level; i++) _octree.SplitAllLeaves();

        }
        /// <summary>
        /// ‚·‚×‚Ä‚Ì‚RŠpŒ`‚ğ“o˜^
        /// </summary>
        public void RegistTriangles()
        {
            for (int i = 0; i < _triangles.Count; i++)
            {
                _octree.RegistTriangle(_triangles[i]);
            }
        }

        public void AddTriangle(Triangle triangle)
        {
            _triangles.Add(triangle);
        }


        // ƒŒƒC‚ğ”ò‚Î‚µ‚ÄAˆê”ÔÅ‰‚É‚Ô‚Â‚©‚é“_‚ğ•Ô‹pB
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor)
        {
            Vector3d dummy1;
            Triangle dummy2;
            IntersectionStatus istat = this.RayTrace(ray, out pos, out nor, out dummy1, out dummy2, CullMode.CullNone);

            return istat;
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan)
        {
            Triangle dummy;
            IntersectionStatus istat = this.RayTrace(ray, out pos, out nor, out tan, out dummy, CullMode.CullNone);

            return istat;
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan, out Triangle tri, CullMode cmode)
        {
            IntersectionStatus istat = _octree.RayTrace(ray, out pos, out nor, out tan, out tri, cmode);

            return istat;
        }
        public OctreeScene()
        {
        }

    }

    //[Serializable()]
    //public class Scene
    //{
    //    // UniformGrid
    //    UniformGrid _ug;

    //    public UniformGrid UGrid
    //    {
    //        get { return _ug; }
    //    }

    //    // ‚RŠpŒ`
    //    private List<Triangle> _triangles = new List<Triangle>();

    //    public Triangle GetTriangle(int idx)
    //    {
    //        return _triangles[idx];
    //    }
    //    public int NumTriangles
    //    {
    //        get { return _triangles.Count; }
    //    }

    //    public Scene()
    //    {
    //    }

    //    public void AddTriangle(Triangle triangle)
    //    {
    //        _triangles.Add(triangle);
    //    }
    //    public void SetupUniformGrid(AABB aabb, int div)
    //    {
    //        _ug = new UniformGrid(aabb, div);
    //    }

    //    /// <summary>
    //    /// ‚·‚×‚Ä‚Ì‚RŠpŒ`‚ğUniformGrid‚É“o˜^
    //    /// </summary>
    //    public void RegistTriangles()
    //    {
    //        for (int i = 0; i < _triangles.Count; i++)
    //        {
    //            //_ug.RegistTriangle(_triangles[i]);
    //            _ug.RegistTriangle(_triangles[i]);
    //        }
    //    }

    //    // ƒŒƒC‚ğ”ò‚Î‚µ‚ÄAˆê”ÔÅ‰‚É‚Ô‚Â‚©‚é“_‚ğ•Ô‹pB
    //    public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor)
    //    {
    //        IntersectionStatus istat = _ug.RayTrace(ray, out pos, out nor);

    //        return istat;
    //    }
    //    public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan)
    //    {
    //        IntersectionStatus istat = _ug.RayTrace(ray, out pos, out nor, out tan);

    //        return istat;
    //    }
    //    public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan, out Triangle tri, CullMode cmode)
    //    {
    //        IntersectionStatus istat = _ug.RayTrace(ray, out pos, out nor, out tan, out tri, cmode);

    //        return istat;
    //    }
    //}
}
