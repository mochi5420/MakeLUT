using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public enum CullMode
    {
        CullNone,
        Clockwise,
        CounterClockwise,
    };

    [Serializable()]
    public class _UniformGrid
    {
        // Cell
        private Cell[] _cells;
        public Cell[] Cell
        {
            get { return _cells; }
        }
        public int NumCells
        {
            get { return _cells.Length; }
        }
        int _div;

        // UniformGrid���̂�Bounding Box
        private AABB _aabb;
        public AABB BoundingBox
        {
            get { return _aabb; }
            set { _aabb = value; }
        }

        public class UniformGridNode
        {
            private AABB _aabb;

            public AABB AABB
            {
                get { return _aabb; }
                set { _aabb = value; }
            }
            private UniformGridNode[] _children;
            private UniformGridNode _parent = null;

            public UniformGridNode Parent
            {
                get { return _parent; }
                set { _parent = value; }
            }

            public UniformGridNode[] Children
            {
                get { return _children; }
                set { _children = value; }
            }


            private UniformGridNode _nextBrother = null;

            public UniformGridNode NextBrother
            {
                get { return _nextBrother; }
                set { _nextBrother = value; }
            }

            public UniformGridNode(UniformGridNode parent) { _parent = parent; }

            public void InsertChild(int NumChild)
            {
                UniformGridNode[] children = new UniformGridNode[NumChild];
                for (int i = 0; i < NumChild; i++) children[i] = new UniformGridNode(this);
                this.Children = children;

                // �Z��o�^
                for (int i = 0; i < NumChild - 1; i++) children[i].NextBrother = children[i + 1];
                children[NumChild - 1].NextBrother = null;

            }

        };
        public class UniformGridOctree
        {
            private int _level;
            public int Level
            {
                get { return _level; }
                set { _level = value; }
            }
            private UniformGridNode _root;

            public UniformGridNode Root
            {
                get { return _root; }
                set { _root = value; }
            }

            /// <summary>
            /// �w��̉ƌn�E�[���̃m�[�h��Ԃ��B
            /// </summary>
            /// <param name="level">�[��</param>
            /// <param name="indices">�ƌn</param>
            /// <returns></returns>
            public UniformGridNode GetNode(int level, int[] indices)
            {
                UniformGridNode node = _root;

                for (int i = 0; i < level; i++)
                {
                    node = node.Children[indices[i]];
                    if (node == null)
                    {
                        throw new Exception("Null Node Exception");
                    }
                }

                return node;
            }


            public void Setup(int level)
            {
                _level = level;
                _root = new UniformGridNode(null);
                UniformGridNode node = _root;
                int current_level = 0;

                while (true)
                {
                    // �q���ǉ�
                    // current_level = n�̂Ƃ���n+1��ڂ̎q�����ǉ�����Ă���B
                    if (current_level == level) break;

                    node.InsertChild(8);

                    // �Z���{��
                    if (node.NextBrother != null)
                    {
                        // �Z�킪����Ƃ�
                        // �������̌Z��ɐݒ�
                        node = node.NextBrother;

                    }
                    else
                    {
                        // �Z�킪���Ȃ��Ƃ��i�����q�̏ꍇ�j
                        // ���j�ɐݒ�
                        node = node.Children[0];
                        current_level++;

                    }
                }
            }
        }

        private UniformGridOctree _octree;

        public _UniformGrid(AABB aabb, int level)
        {
            _aabb = aabb;

            Vector3d MinVec = aabb.Center - aabb.HalfVector;
            Vector3d MaxVec = aabb.Center + aabb.HalfVector;

            _div = level * level * level;

            // octree�\�z
            //_octree = new UniformGridOctree();
            //_octree.Setup(level);
            // �����܂�

            Setup(MinVec.X, MaxVec.X, MinVec.Y, MaxVec.Y, MinVec.Z, MaxVec.Z, _div);
            //SetupNew(MinVec.X, MaxVec.X, MinVec.Y, MaxVec.Y, MinVec.Z, MaxVec.Z, _div);
        }
        //public UniformGrid(AABB aabb, int div)
        //{
        //    _aabb = aabb;

        //    Vector3d MinVec = aabb.Center - aabb.HalfVector;
        //    Vector3d MaxVec = aabb.Center + aabb.HalfVector;

        //    _div = div;

        //    Setup(MinVec.X, MaxVec.X, MinVec.Y, MaxVec.Y, MinVec.Z, MaxVec.Z, div);
        //}
        public int CellId(int xid, int yid, int zid)
        {
            if (xid < 0 || xid >= _div) return -1;
            if (yid < 0 || yid >= _div) return -1;
            if (zid < 0 || zid >= _div) return -1;

            return _div * _div * xid + _div * yid + zid;
        }
        // ���alen��Cell��Ԃ�
        public Cell[] GetCells(Vector3d pos, double len)
        {
            Vector3d u = 2.0 / (double)_div * _aabb.HalfVector;

            int CellXNum = (int)Math.Ceiling(len / u.X);
            int CellYNum = (int)Math.Ceiling(len / u.Y);
            int CellZNum = (int)Math.Ceiling(len / u.Z);

            List<Cell> CellList = new List<Cell>();
            int id_x, id_y, id_z;
            GetCellID(pos, out id_x, out id_y, out id_z);  // ����

            for (int i = -CellXNum; i <= CellXNum; i++)
            {
                for (int j = -CellYNum; j <= CellYNum; j++)
                {
                    for (int k = -CellZNum; k <= CellZNum; k++)
                    {
                        int id = CellId(id_x + i, id_y + j, id_z + k);
                        if (id != -1) CellList.Add(Cell[id]);
                    }
                }
            }

            Cell[] cells = new Cell[CellList.Count];
            for (int i = 0; i < CellList.Count; i++) cells[i] = CellList[i];

            return cells;
        }
        public Cell[][] GetCellsLevel(Vector3d pos)
        {
            Vector3d u = 2.0 / (double)_div * _aabb.HalfVector;

            List<Cell> CellList = new List<Cell>();

            List<int>[] IDList = new List<int>[_div];

            int id_x, id_y, id_z;
            GetCellID(pos, out id_x, out id_y, out id_z);  // ����

            id_x = 10;
            id_y = 0;
            id_z = 0;

            for (int level = 0; level < _div; level++)
            {
                IDList[level] = new List<int>();
                for (int i = 0; i < 6; i++)
                {
                    for (int j = -level; j <= level; j++)
                    {
                        for (int k = -level; k <= level; k++)
                        {
                            int x = 0, y = 0, z = 0;

                            if (i == 0) { x = id_x + level; y = id_y + j; z = id_z + k; }
                            else if (i == 1) { x = id_x - level; y = id_y + j; z = id_z + k; }
                            else if (i == 2) { x = id_x + k; y = id_y + level; z = id_z + i; }
                            else if (i == 3) { x = id_x + k; y = id_y - level; z = id_z + i; }
                            else if (i == 4) { x = id_x + i; y = id_y + j; z = id_z + level; }
                            else if (i == 5) { x = id_x + i; y = id_y + j; z = id_z - level; }

                            int id = CellId(x, y, z);
                            if (id != -1) IDList[level].Add(id);
                        }
                    }
                }
            }

            int MaxLevel = 0;
            for (int i = 0; i < IDList.Length; i++)
            {
                if (IDList[i].Count != 0) MaxLevel = i;
            }

            Cell[][] cells = new Cell[MaxLevel][];
            for (int i = 0; i < MaxLevel; i++)
            {
                cells[i] = new Cell[IDList[i].Count];
                for (int j = 0; j < IDList[i].Count; j++)
                {
                    int id = IDList[i][j];
                    cells[i][j] = Cell[id];
                }
            }

            return cells;
        }

        // pos���܂�Cell��ID��Ԃ�
        public int GetCellID(Vector3d pos)
        {
            int xid, yid, zid;
            return GetCellID(pos, out xid, out yid, out zid);
        }
        // pos���܂�Cell��ID��Ԃ�
        public int GetCellID(Vector3d pos, out int xid, out int yid, out int zid)
        {
            Vector3d u = 2.0 / (double)_div * _aabb.HalfVector;
            Vector3d cpos = pos - _aabb.Center + _aabb.HalfVector;

            int x_idx = (int)Math.Floor(cpos.X / u.X);
            int y_idx = (int)Math.Floor(cpos.Y / u.Y);
            int z_idx = (int)Math.Floor(cpos.Z / u.Z);

            if (x_idx < 0) x_idx = 0;
            else if (x_idx >= _div) x_idx = _div - 1;
            if (y_idx < 0) y_idx = 0;
            else if (y_idx >= _div) y_idx = _div - 1;
            if (z_idx < 0) z_idx = 0;
            else if (z_idx >= _div) z_idx = _div - 1;

            xid = x_idx;
            yid = y_idx;
            zid = z_idx;
            return _div * _div * x_idx + _div * y_idx + z_idx;
        }
 

        public void Setup(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int div)
        {
            // octree�\�z
            _octree = new UniformGridOctree();
            _octree.Setup(div);
            // �����܂�

            double ux = (x_max - x_min) / (double)div;
            double uy = (y_max - y_min) / (double)div;
            double uz = (z_max - z_min) / (double)div;

            Vector3d boxhalf = 0.5 * new Vector3d(ux, uy, uz);

            _cells = new Cell[div * div * div];

            for (int i = 0; i < div; i++)
            {
                for (int j = 0; j < div; j++)
                {
                    for (int k = 0; k < div; k++)
                    {
                        Vector3d center = new Vector3d((double)i * ux + boxhalf.X + x_min, (double)j * uy + boxhalf.Y + y_min, (double)k * uz + boxhalf.Z + z_min);

                        _cells[div * div * i + div * j + k] = new Cell(new AABB(center, boxhalf));

                    }
                }
            }
        }

        public void SetupNew(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int div)
        {

            int level = _octree.Level;

            UniformGridNode node = _octree.Root;

            for (int cl = 0; cl < level; cl++)
            {
                int d = 1;
                for (int i = 0; i < cl; i++) d *= 2;    // 2�ׂ̂���

                double ux = (x_max - x_min) / (double)d;
                double uy = (y_max - y_min) / (double)d;
                double uz = (z_max - z_min) / (double)d;

                Vector3d boxhalf = 0.5 * new Vector3d(ux, uy, uz);
                
                //_cells = new Cell[div * div * div];

                for (int i = 0; i < d; i++)
                {
                    for (int j = 0; j < d; j++)
                    {
                        for (int k = 0; k < d; k++)
                        {

                            Vector3d center = new Vector3d((double)i * ux + boxhalf.X + x_min, (double)j * uy + boxhalf.Y + y_min, (double)k * uz + boxhalf.Z + z_min);
                            AABB aabb = new AABB(center, boxhalf);

                            node.AABB = aabb;

                            if (node.NextBrother != null)
                            {
                                node = node.NextBrother;
                            }
                            else
                            {
                                node = node.Children[0];
                            }

                            // �����cl�ł킩��B
                            // �ƌn�𒲂ׂȂ���B

                            //_cells[div * div * i + div * j + k] = new Cell(new AABB(center, boxhalf));

                        }
                    }
                }

                //node = node.Children[0];
            }
        }

        // Triangle��AABB�̓����Ɋ܂܂��S����CELL�ƌ���������s���B
        public void RegistTriangle(Triangle tri)
        {
            Vector3d MinVec = new Vector3d(
                Math.Min(Math.Min(tri.V0.X, tri.V1.X), tri.V2.X),
                Math.Min(Math.Min(tri.V0.Y, tri.V1.Y), tri.V2.Y),
                Math.Min(Math.Min(tri.V0.Z, tri.V1.Z), tri.V2.Z));
            Vector3d MaxVec = new Vector3d(
                Math.Max(Math.Max(tri.V0.X, tri.V1.X), tri.V2.X),
                Math.Max(Math.Max(tri.V0.Y, tri.V1.Y), tri.V2.Y),
                Math.Max(Math.Max(tri.V0.Z, tri.V1.Z), tri.V2.Z));

            // Triangle��AABB
            AABB TriAABB = new AABB(0.5 * (MaxVec + MinVec), 0.5 * (MaxVec - MinVec));

            int min_xid, min_yid,min_zid;
            int max_xid, max_yid,max_zid;
            GetCellID(TriAABB.Center - TriAABB.HalfVector, out min_xid, out min_yid, out min_zid);
            GetCellID(TriAABB.Center + TriAABB.HalfVector, out max_xid, out max_yid, out max_zid);

            for (int i = min_xid; i <= max_xid; i++)
            {
                for (int j = min_yid; j <= max_yid; j++)
                {
                    for (int k = min_zid; k <= max_zid; k++)
                    {
                        int cell_id = CellId(i, j, k);

                        if (Intersections.triBoxOverlap(_cells[cell_id].Aabb, tri) == OverlapStatus.Disjoint) continue;
                        else
                        {
                            // 3�p�`��ǉ�
                            _cells[cell_id].TryAddTriangle(tri);
                        }

                    }
                }
            }
        }
        // �SCELL�ƌ���������s���B
        public void RegistTriangleOld(Triangle tri)
        {
            Vector3d MinVec = new Vector3d(
                Math.Min(Math.Min(tri.V0.X, tri.V1.X), tri.V2.X),
                Math.Min(Math.Min(tri.V0.Y, tri.V1.Y), tri.V2.Y),
                Math.Min(Math.Min(tri.V0.Z, tri.V1.Z), tri.V2.Z));
            Vector3d MaxVec = new Vector3d(
                Math.Max(Math.Max(tri.V0.X, tri.V1.X), tri.V2.X),
                Math.Max(Math.Max(tri.V0.Y, tri.V1.Y), tri.V2.Y),
                Math.Max(Math.Max(tri.V0.Z, tri.V1.Z), tri.V2.Z));

            // Triangle��AABB
            AABB TriAABB = new AABB(0.5 * (MaxVec + MinVec), 0.5 * (MaxVec - MinVec));


            for (int i = 0; i < _cells.Length; i++)
            {

                // ��������AABB���������Ȃ��Ȃ玟��CELL�ɐi�ށB
                //if (Intersections.AABBOverlap(TriAABB, _cells[i].Aabb) == OverlapStatus.Disjoint) continue;
                //else
                //{
                //    // 3�p�`��ǉ�
                //    _cells[i].TryAddTriangle(tri);
                //}

                if (Intersections.triBoxOverlap(_cells[i].Aabb, tri) == OverlapStatus.Disjoint) continue;
                else
                {
                    // 3�p�`��ǉ�
                    _cells[i].TryAddTriangle(tri);
                }
            }
        }
        //
        public int[] RayOverlapCells(Ray ray)
        {
            List<int> list = new List<int>();
            
            // ���ׂĂ�Cell��Ray�̃I�[�o�[���b�v����
            for (int i = 0; i < _cells.Length; i++)
            {
                if (_cells[i].RayAABBOverlap(ray) == OverlapStatus.Overlap) list.Add(i);
            }

            // �f�[�^���ڂ�
            int[] indices = new int[list.Count];
            for (int i = 0; i < list.Count; i++)
            {
                indices[i] = (int)list[i];
            }

            return indices;            
        }
        public int[] RayOverlapCellsNew(Ray ray)
        {
            List<int> list = new List<int>();

            // ���ׂĂ�Cell��Ray�̃I�[�o�[���b�v����
            for (int i = 0; i < _cells.Length; i++)
            {
                if (_cells[i].RayAABBOverlap(ray) == OverlapStatus.Overlap) list.Add(i);
            }

            // �f�[�^���ڂ�
            int[] indices = new int[list.Count];
            for (int i = 0; i < list.Count; i++)
            {
                indices[i] = (int)list[i];
            }

            return indices;
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor)
        {
            Vector3d tan;
            IntersectionStatus istat = RayTrace(ray, out  pos, out  nor, out tan);

            return istat;
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan)
        {
            Triangle dummy;
            return RayTrace(ray, out pos, out nor, out tan, out dummy, CullMode.CullNone);
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan, CullMode cmode)
        {
            Triangle dummy;
            return RayTrace(ray, out pos, out nor, out tan, out dummy, cmode);
        }
        public IntersectionStatus RayTrace(Ray ray, out Vector3d pos, out Vector3d nor, out Vector3d tan, out Triangle triangle, CullMode cmode)
        {
            IntersectionStatus istat = IntersectionStatus.Reject;
            triangle = null;

            // �I�[�o�[���b�v����Cell��������
            int[] list = RayOverlapCells(ray);

            // ���̃Z�������R�p�`�ƏՓ˔���
            double min_t = double.MaxValue;
            pos = new Vector3d();
            nor = new Vector3d();
            tan = new Vector3d();
            for (int i = 0; i < list.Length; i++)
            {
                for (int j = 0; j < _cells[list[i]].NumTriangle; j++)
                {
                    Triangle tri = _cells[list[i]].GetTriangle(j);
                    double u, v, t;
                    Vector3d p;

                    // �Փ˔���
                    if (Intersections.RayTriangleIntersect(ray, tri, out u, out v, out t, out p) == IntersectionStatus.Intersect)
                    {
                        // �����Triangle�ɂԂ�������C�ɂ��Ȃ��B
                        if (t < 0.0) continue;

                        // �������̖ʂ͋C�ɂ��Ȃ��B

                        if (cmode == CullMode.Clockwise)
                        {
                            if (Vector3d.Dot(ray.Direction, tri.FaceNormal) > 0.0) continue;
                        }
                        else if (cmode == CullMode.CounterClockwise)
                        {
                            if (Vector3d.Dot(ray.Direction, tri.FaceNormal) < 0.0) continue;
                        }

                        // ����O��Triangle�ɂԂ����Ă�����B
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

    // �O���b�h�̂P���Z��
    [Serializable()]
    public class Cell
    {
        // Cell�͈̔͂��߂�AABB
        private AABB _aabb;

        public AABB Aabb
        {
            get { return _aabb; }
            set { _aabb = value; }
        }

        // Cell�ɏ�������O�p�`
        private List<Triangle> _triangles = new List<Triangle>();

        // �o�^���ꂽ3�p�`�̌���Ԃ�
        public int NumTriangle
        {
            get { return _triangles.Count; }
        }

        public Triangle GetTriangle(int idx)
        {
            return (Triangle)_triangles[idx];
        }


        // ��X�ɂ̓t�H�g�����i�[�ł������������ȁB

        // �R���X�g���N�^
        public Cell(AABB aabb)
        {
            _aabb = aabb;
        }


        // 3�p�`���I�[�o�[���b�v���Ă�����ǉ�����B
        // �ǉ����ꂽ��true, ����Ȃ�������false��Ԃ��B
        public bool TryAddTriangle(Triangle tri)
        {
            if (Intersections.triBoxOverlap(this._aabb, tri) == OverlapStatus.Overlap)
            {
                // 3�p�`���̂��̂������Ă�̂Ń��������������邩���B
                _triangles.Add(tri);
                return true;
            }
            else
            {
                return false;
            }
        }

        // ����Ray�i�����j�Ƃ��̃Z����AABB�Ƃ̃I�[�o�[���b�v����
        public OverlapStatus RayAABBOverlap(Ray ray)
        {
            return Intersections.RayAABBOverlap(ray, this._aabb);
        }

        //// ����Ray�i�����j�ƁA���̃Z���ɓo�^����Ă���R�p�`�̏Փ˂𔻒�B
        //public IntersectionStatus RayTriangleIntersection(Ray ray, out Vector3d pos)
        //{
        //    IntersectionStatus istat = IntersectionStatus.Reject;
        //    pos = new Vector3d();

        //    double min_t = double.MaxValue;

        //    // ��Ԏ�O�łԂ���Triangle��T���B
        //    for (int i = 0; i < _triangles.Count; i++)
        //    {
        //        double u, v, t;
        //        Vector3d p;

        //        // �Փ˂���ꏊ�����߂�B
        //        if (Intersections.RayTriangleIntersect(ray, (Triangle)_triangles[i], out u, out v, out t, out p) == IntersectionStatus.Intersect)
        //        {
        //            // ���łԂ����Ă��疳���B
        //            if (t < 0.0) continue;

        //            // ����O�łԂ����Ă�����X�V�B
        //            if (t < min_t)
        //            {
        //                min_t = t;
        //                pos = p;

        //                istat = IntersectionStatus.Intersect;
        //            }
        //        }
        //    }

        //    return istat;
        //}
    }

}
