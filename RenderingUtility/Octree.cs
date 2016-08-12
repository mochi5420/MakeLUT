using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using Utility;

namespace RenderingUtility
{
    public class Octree
    {
        private OctreeNode _root;

        public OctreeNode Root
        {
            get { return _root; }
        }

        public Octree()
            : this(1, -1, 1, -1, 1, -1)
        {
        }
        public Octree(AABB aabb)
            : this(aabb.Right, aabb.Left, aabb.Front, aabb.Back, aabb.Top, aabb.Bottom)
        {
        }

        public Octree(double xMax, double xMin, double yMax, double yMin, double zMax, double zMin)
        {
            _root = new OctreeNode(xMax, xMin, yMax, yMin, zMax, zMin, "Root");
        }

        public static void PreOrder(OctreeNode node)
        {
            OctreeNode thisNode = node;

            // 処理
            System.Diagnostics.Debug.WriteLine(node.UserObject.ToString());
            //

            if (node.hasChildren())
            {
                for (int i = 0; i < node.Branch.Length; i++)
                {
                    PreOrder(node.Branch[i]);
                }
            }
        }


        public void SplitOne(int level)
        {
            OctreeNode node = _root;
            for (int i = 0; i < level; i++)
            {
                node = node.Branch[0];
            }
            node.Split();
        }   
    

        public void SplitAllLeaves()
        {
            Split(_root);
        }

        private static void Split(OctreeNode node)
        {
            OctreeNode thisNode = node;

            if (node.hasChildren())
            {
                for (int i = 0; i < node.Branch.Length; i++)
                {
                    Split(node.Branch[i]);
                }
            }

            // 帰りがけに処理（子がいないノードならSplit）
            if (!node.hasChildren()) node.Split();
            //

        }
        // 葉をすべて取得
        public OctreeNode[] GetAllLeaves()
        {
            List<OctreeNode> list = new List<OctreeNode>();
            GetLeaf(_root, list);

            OctreeNode[] nodes = new OctreeNode[list.Count];
            for (int i = 0; i < nodes.Length; i++) nodes[i] = list[i];

            return nodes;
        }
        private void GetLeaf(OctreeNode node, List<OctreeNode> list)
        {
            if (node.hasChildren())
            {
                // 子供を調査
                for (int i = 0; i < node.Branch.Length; i++)
                {
                    GetLeaf(node.Branch[i], list);
                }
            }
            else
            {
                list.Add(node);
            }
        }
        public OctreeNode[] GetOverlapLeaves(Vector3d center, double radius)
        {
            List<OctreeNode> list = new List<OctreeNode>();
            GetOverlapLeaves(_root, center, radius, list);

            OctreeNode[] nodes = new OctreeNode[list.Count];
            for (int i = 0; i < nodes.Length; i++) nodes[i] = list[i];

            return nodes;
        }
        public OctreeNode[] GetTriangleOverlapLeaves(Triangle triangle)
        {
            List<OctreeNode> list = new List<OctreeNode>();
            GetTriangleOverlapLeaves(_root, triangle, list);

            OctreeNode[] nodes = new OctreeNode[list.Count];
            for (int i = 0; i < nodes.Length; i++) nodes[i] = list[i];

            return nodes;
        }
        private void GetPointOverlapLeaf(OctreeNode node, Vector3d pos, List<OctreeNode> list)
        {
            // 自分とPointのオーバーラップをテスト
            OverlapStatus ostat = Intersections.PointAABBOverlap(pos, node.Bounds);

            // debug
            //list.Add(node);

            if (ostat == OverlapStatus.Disjoint)
            {
                // もしオーバーラップしなかったら（子供がいてもいなくても）
                // 子孫はオーバーラップしないから打ち切り。
                return;
            }
            else if (node.hasChildren())
            {
                // もしオーバーラップしてて子供がいれば、子供を調査
                for (int i = 0; i < node.Branch.Length; i++)
                {
                    GetPointOverlapLeaf(node.Branch[i], pos, list);
                }

            }
            else
            {
                // もしオーバーラップしてて子供がいなければリストに追加
                list.Add(node);
                return;
            }
        }
        
        private void GetOverlapLeaves(OctreeNode node, Vector3d center, double radius, List<OctreeNode> list)
        {
            // 自分とSphereのオーバーラップをテスト
            OverlapStatus ostat = Intersections.SphereAABBOverlap(center, radius, node.Bounds);


            if (ostat == OverlapStatus.Disjoint)
            {
                // もしオーバーラップしなかったら（子供がいてもいなくても）
                // 子孫はオーバーラップしないから打ち切り。
                return;
            }
            else if (!node.hasChildren())
            {
                // もしオーバーラップしてて子供がいなければリストに追加
                list.Add(node);
            }
            else
            {
                // もしオーバーラップしてて子供がいれば、子供を調査

                for (int i = 0; i < node.Branch.Length; i++)
                {
                    GetOverlapLeaves(node.Branch[i], center, radius, list);
                }

            }
        }



        private void GetTriangleOverlapLeaves(OctreeNode node, Triangle triangles, List<OctreeNode> list)
        {
            // 自分とTriangleのオーバーラップをテスト
            OverlapStatus ostat = Intersections.triBoxOverlap(node.Bounds, triangles);

            if (ostat == OverlapStatus.Disjoint)
            {
                // もしオーバーラップしなかったら（子供がいてもいなくても）
                // 子孫はオーバーラップしないから打ち切り。
                return;
            }
            else if (!node.hasChildren())
            {
                // もしオーバーラップしてて子供がいなければリストに追加
                list.Add(node);
            }
            else
            {
                // もしオーバーラップしてて子供がいれば、子供を調査

                for (int i = 0; i < node.Branch.Length; i++)
                {
                    GetTriangleOverlapLeaves(node.Branch[i], triangles, list);
                }
            }
        }
        private void GetRayOverlapLeaves(OctreeNode node, Ray ray, List<OctreeNode> list)
        {
            // 自分とTriangleのオーバーラップをテスト
            OverlapStatus ostat = Intersections.RayAABBOverlap(ray, node.Bounds);

            if (ostat == OverlapStatus.Disjoint)
            {
                // もしオーバーラップしなかったら（子供がいてもいなくても）
                // 子孫はオーバーラップしないから打ち切り。
                return;
            }
            else if (!node.hasChildren())
            {
                // もしオーバーラップしてて子供がいなければリストに追加
                list.Add(node);
            }
            else
            {
                // もしオーバーラップしてて子供がいれば、子供を調査

                for (int i = 0; i < node.Branch.Length; i++)
                {
                    GetRayOverlapLeaves(node.Branch[i], ray, list);
                }
            }
        }

        public void RegistTriangle(Triangle triangle)
        {
            List<OctreeNode> list = new List<OctreeNode>();

            // Triangleと重なるノードを探索
            GetTriangleOverlapLeaves(_root, triangle, list);

            // 重なったノードにTriangleを登録
            for (int i = 0; i < list.Count; i++) list[i].AddTriangle(triangle);

        }

        public void RegistTriangle2(Triangle triangle)
        {
            List<OctreeNode> list = new List<OctreeNode>();

            // Triangleと重なるノードを探索
            GetTriangleOverlapLeaves(_root, triangle, list);

            // 重なったノードにTriangleを登録
            for (int i = 0; i < list.Count; i++)
            {
                list[i].AddTriangle(triangle);
            }
        }

        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan, out Triangle triangle, CullMode cmode)
        {
            IntersectionStatus istat = IntersectionStatus.Reject;
            triangle = null;

            
            // オーバーラップするCellを見つける
            List<OctreeNode> list = new List<OctreeNode>();
            GetRayOverlapLeaves(_root, ray, list);

            {
            //Intersections.RayAABBIntersect(

            //for(int i=0; i<list.Count; i++) list[i].ValSort = Vector3d.LengthSq(ray.Origin-

            //list.Sort();
            }
            
            // そのセルがもつ３角形と衝突判定
            double min_t = double.MaxValue;
            pos = new Vector3d();
            nor = new Vector3d();
            tan = new Vector3d();
            for (int i = 0; i < list.Count; i++)
            {
                for (int j = 0; j < list[i].Triangles.Count; j++)
                {
                    Triangle tri = list[i].Triangles[j];
                    double u, v, t;
                    Vector3d p;

                    // 衝突判定
                    if (Intersections.RayTriangleIntersect(ray, tri, out u, out v, out t, out p) == IntersectionStatus.Intersect)
                    {
                        // 後方のTriangleにぶつかったら気にしない。
                        if (t < 0.0) continue;

                        // 裏向きの面は気にしない。

                        if (cmode == CullMode.Clockwise)
                        {
                            if (Vector3d.Dot(ray.Direction, tri.FaceNormal) > 0.0) continue;
                        }
                        else if (cmode == CullMode.CounterClockwise)
                        {
                            if (Vector3d.Dot(ray.Direction, tri.FaceNormal) < 0.0) continue;
                        }

                        // より手前のTriangleにぶつかっていたら。
                        if (t < min_t)
                        {
                            min_t = t;
                            pos = p;
                            //nor = tri.Normal;
                            nor = Vector3d.Normalize(tri.Normal[0] + u * (tri.Normal[1] - tri.Normal[0]) + v * (tri.Normal[2] - tri.Normal[0]));
                            tan = Vector3d.Normalize(tri.Tangent[0] + u * (tri.Tangent[1] - tri.Tangent[0]) + v * (tri.Tangent[2] - tri.Tangent[0]));
                            triangle = tri;
                            istat = IntersectionStatus.Intersect;
                        }
                    }
                }
            }
            
            
            return istat;
        }
         
    }

    public class OctreeNode : IComparable
    {
        private AABB _bounds;
        private OctreeNode[] _branch;
        private Object _object;

        private List<Triangle> _triangles = new List<Triangle>();


        private double _ValSort;

        public double ValSort
        {
            get { return _ValSort; }
            set { _ValSort = value; }
        }

        public List<Triangle> Triangles
        {
            get { return _triangles; }
            set { _triangles = value; }
        }

        public OctreeNode(double xMax, double xMin, double yMax, double yMin, double zMax, double zMin)
            : this(xMax, xMin, yMax, yMin, zMax, zMin, null)
        {
        }
        public OctreeNode(double xMax, double xMin, double yMax, double yMin, double zMax, double zMin, Object obj)
        {
            _bounds = new AABB(xMax, xMin, yMax, yMin, zMax, zMin);
            _object = obj;
        }


        // 比較
        public int CompareTo(object other)
        {
            return this._ValSort.CompareTo(((OctreeNode)other)._ValSort);
        }

        public OctreeNode[] Branch
        {
            get { return _branch; }
        }
        //public OctreeNode Parent
        //{
        //    get { return _parent; }
        //}
        public Object UserObject
        {
            get { return _object; }
        }
        public AABB Bounds
        {
            get { return _bounds; }
        }

        public void SetUserObject(string text)
        {
            _object = text;
        }

        /// <summary>Return true if the node has branch. </summary>
        public bool hasChildren()
        {
            if (_branch != null)
                return true;
            else
                return false;
        }



        public void Split()
        {
            float nsHalf = (float)(_bounds.Top - (_bounds.Top - _bounds.Bottom) * 0.5);
            float ewHalf = (float)(_bounds.Right - (_bounds.Right - _bounds.Left) * 0.5);
            float fbHalf = (float)(_bounds.Front - (_bounds.Front - _bounds.Back) * 0.5);

            _branch = new OctreeNode[8];

            string parentName = _object.ToString();

            _branch[0] = new OctreeNode(ewHalf, _bounds.Left, _bounds.Front, fbHalf, _bounds.Top, nsHalf, parentName + "_0"); //left-front-top
            _branch[1] = new OctreeNode(_bounds.Right, ewHalf, _bounds.Front, fbHalf, _bounds.Top, nsHalf, parentName + "_1");
            _branch[2] = new OctreeNode(ewHalf, _bounds.Left, _bounds.Front, fbHalf, nsHalf, _bounds.Bottom, parentName + "_2");
            _branch[3] = new OctreeNode(_bounds.Right, ewHalf, _bounds.Front, fbHalf, nsHalf, _bounds.Bottom, parentName + "_3");

            _branch[4] = new OctreeNode(ewHalf, _bounds.Left, fbHalf, _bounds.Back, _bounds.Top, nsHalf, parentName + "_4");//left-back-top
            _branch[5] = new OctreeNode(_bounds.Right, ewHalf, fbHalf, _bounds.Back, _bounds.Top, nsHalf, parentName + "_5");
            _branch[6] = new OctreeNode(ewHalf, _bounds.Left, fbHalf, _bounds.Back, nsHalf, _bounds.Bottom, parentName + "_6");
            _branch[7] = new OctreeNode(_bounds.Right, ewHalf, fbHalf, _bounds.Back, nsHalf, _bounds.Bottom, parentName + "_7");
        }
        public void AddTriangle(Triangle tri)
        {
            _triangles.Add(tri);
        }

        
    }

}
