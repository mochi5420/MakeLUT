using System;
using System.Collections.Generic;
using System.Text;
using MathUtil;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace MeshUtil
{
    public class WavefrontObj
    {
        public enum Format
        {
            Ascii,
            OriginalBinary
        };

        public class Facet
        {
            private int[] _posIdx;
            private int[] _norIdx;
            private int[] _tcIdx;
            private int[] _tanIdx;

            public int[] PosIdx { get { return _posIdx; } set { _posIdx = value; } }
            public int[] NorIdx { get { return _norIdx; } set { _norIdx = value; } }
            public int[] TCIdx { get { return _tcIdx; } set { _tcIdx = value; } }
            public int[] TanIdx { get { return _tanIdx; } set { _tanIdx = value; } }

            private int _materialId;

            public int MaterialId
            {
                get { return _materialId; }
                set { _materialId = value; }
            }

            private int _groudId = 0;
            public int GroupId { get { return _groudId; } set { _groudId = value; } }

            /// <summary>
            /// 何角形かを返す
            /// </summary>
            public int NumEdge { get { return _posIdx.Length; } }


            public bool IsVertex
            {
                get { return _posIdx != null ? true : false; }

            }

            public bool IsNormal
            {
                get
                {
                    if (_norIdx == null) return false;
                    foreach (int idx in _norIdx)
                    {
                        if (idx == -1) return false;
                    }
                    return true;
                }
            }

            public bool IsTexCoord
            {
                get
                {
                    if (_tcIdx == null) return false;
                    foreach (int idx in _tcIdx)
                    {
                        if (idx == -1) return false;
                    }
                    return true;
                }
            }
            public bool IsTangent
            {
                get
                {
                    if (_tanIdx == null) return false;
                    foreach (int idx in _tanIdx)
                    {
                        if (idx == -1) return false;
                    }
                    return true;
                }
            }
            public Facet[] Triangulate()
            {
                Facet[] facets = new Facet[this.NumEdge - 2];

                for (int i = 0; i < this.NumEdge - 2; i++)
                {
                    int[] pidx = null;
                    int[] nidx = null;
                    int[] tcidx = null;

                    if (IsVertex)
                    {
                        pidx = new int[3];
                        pidx[0] = this.PosIdx[0];
                        pidx[1] = this.PosIdx[i + 1];
                        pidx[2] = this.PosIdx[i + 2];
                    }

                    if (IsNormal)
                    {
                        nidx = new int[3];
                        nidx[0] = this.NorIdx[0];
                        nidx[1] = this.NorIdx[i + 1];
                        nidx[2] = this.NorIdx[i + 2];
                    }

                    if (IsTexCoord)
                    {
                        tcidx = new int[3];
                        tcidx[0] = this.TCIdx[0];
                        tcidx[1] = this.TCIdx[i + 1];
                        tcidx[2] = this.TCIdx[i + 2];
                    }


                    facets[i] = new Facet(pidx, nidx, tcidx);

                }

                return facets;
            }
            public void SetupTangentIdx(int n)
            {
                _tanIdx = new int[n];
            }

            /// <summary>
            /// コンストラクタ
            /// </summary>
            /// <param name="n">n角形</param>
            public Facet(int n)
            {
                _posIdx = new int[n];
                _norIdx = new int[n];
                _tcIdx = new int[n];
            }



            public Facet(int[] pidx, int[] nidx, int[] tcidx)
            {
                if (pidx != null)
                {
                    _posIdx = new int[pidx.Length];
                    for (int i = 0; i < _posIdx.Length; i++) _posIdx[i] = pidx[i];
                }

                if (nidx != null)
                {
                    _norIdx = new int[nidx.Length];
                    for (int i = 0; i < _norIdx.Length; i++) _norIdx[i] = nidx[i];
                }

                if (tcidx != null)
                {
                    _tcIdx = new int[tcidx.Length];
                    for (int i = 0; i < _tcIdx.Length; i++) _tcIdx[i] = tcidx[i];
                }

            }

            public Facet(string[] token)
            {
                // 角形
                int n = token.Length - 1;

                int[] vidx = new int[n];
                int[] vnidx = new int[n];
                int[] vtidx = new int[n];

                bool exist_v = true;
                bool exist_vn = true;
                bool exist_vt = true;

                for (int i = 0; i < n; i++)
                {
                    string[] p = token[i + 1].Split(new char[] { '/' }, StringSplitOptions.None);

                    if (p.Length == 3)  /* v/vt/vn or v//vn */
                    {
                        if (!int.TryParse(p[0], out vidx[i])) { vidx[i] = -1; exist_v = false; }
                        if (!int.TryParse(p[1], out vtidx[i])) { vtidx[i] = -1; exist_vt = false; }
                        if (!int.TryParse(p[2], out vnidx[i])) { vnidx[i] = -1; exist_vn = false; }
                    }
                    else if (p.Length == 2) /* v/vt */
                    {
                        if (!int.TryParse(p[0], out vidx[i])) { vidx[i] = -1; exist_v = false; }
                        if (!int.TryParse(p[1], out vtidx[i])) { vtidx[i] = -1; exist_vt = false; }
                    }
                    else if (p.Length == 1)  /* v */
                    {
                        if (!int.TryParse(p[0], out vidx[i])) { vidx[i] = -1; exist_v = false; }
                    }
                }

                if (exist_v)
                {
                    _posIdx = new int[n];
                    for (int i = 0; i < n; i++) _posIdx[i] = vidx[i] - 1;
                }

                if (exist_vn)
                {
                    _norIdx = new int[n];
                    for (int i = 0; i < n; i++) _norIdx[i] = vnidx[i] - 1;
                }

                if (exist_vt)
                {
                    _tcIdx = new int[n];
                    for (int i = 0; i < n; i++) _tcIdx[i] = vtidx[i] - 1;
                }
            }
        }

        private Vector3d[] _position;
        private Vector3d[] _normal;
        private Vector3d[] _tangent;
        private Vector2d[] _texCoord;
        private Facet[] _facets;

        private MaterialLibrary _mtllib;
        public int NumMaterials { get { return _mtllib.NumMaterials; } }
        public WavefrontMaterial[] Materials { get { return _mtllib.Materials; } }
        public bool IsMaterialLibrary { get { return _mtllib != null ? true : false; } }

        public int NumMaps { get { return _mtllib.NumMaps; } }
        public string MapPath(int mapId) { return _dir + "\\" + _mtllib.MapNames[mapId]; }
        public int MapId(int matId) { return _mtllib.Materials[matId].MapKdId; }
        public int MapKdId(int matId) { return _mtllib.Materials[matId].MapKdId; }
        public int MapBumpId(int matId) { return _mtllib.Materials[matId].MapBumpId; }

        // マテリアルごとのFacetの通し番号
        private int[][] _facetIdPermat;

        private string[] _groupNames;

        public int GetNumMaterialFacet(int matid) { return _facetIdPermat[matid].Length; }
        public int GetMaterialFacet(int matid, int fidx) { return _facetIdPermat[matid][fidx]; }
        public int GetMaterialIdFromTechnique(string tech)
        {
            for (int i = 0; i < _mtllib.Materials.Length; i++) if (tech == _mtllib.Materials[i].Technique) return i;
            return -1;
        }
        public int GetMaterialIdFromName(string mat)
        {
            for (int i = 0; i < _mtllib.Materials.Length; i++) if (mat == _mtllib.Materials[i].Name) return i;
            return -1;
        }



        // マテリアル操作
        public void SetKdByName(string name, Vector3d kd)
        {
            _mtllib.SetKdByName(name, kd);
        }
        public Vector3d GetKdByName(string name)
        {
            return _mtllib.GetKdByName(name);
        }
        public string GetMaterialName(int id)
        {
            return _mtllib.Materials[id].Name;
        }

        // GroupId
        public int GetGroupId(int fidx)
        {
            return _facets[fidx].GroupId;
        }
        public string GetGroupName(int fidx)
        {
            return _groupNames[GetGroupId(fidx)];
        }


        private string _dir;

        public void TransparentMaterialSort()
        {
            _mtllib.TransparentMaterialSort();
        }

        public WavefrontObj()
        {
        }
        public WavefrontObj(string path)
            : this(path, null, Format.Ascii)
        {
        }
        public WavefrontObj(string objpath, string mtlpath)
            : this(objpath, mtlpath, Format.Ascii)
        {
        }
        public WavefrontObj(string objpath, string mtlpath, Format format)
        {
            _dir = Path.GetDirectoryName(objpath);

            // マテリアルが指定されていないとき用に標準のマテリアルを用意
            //_mtllib = new MaterialLibrary();

            if (format == Format.Ascii)
            {
                if (mtlpath != null)
                {
                    _mtllib = new MaterialLibrary(mtlpath);
                    if (!Load(objpath, true)) throw new Exception("ファイルの読み込みに失敗した。");
                }
                else
                {
                    if (!Load(objpath, false)) throw new Exception("ファイルの読み込みに失敗した。");
                }
            }
            else
            {
                _mtllib = new MaterialLibrary(mtlpath);
                if (!LoadBinary(objpath)) throw new Exception("ファイルの読み込みに失敗した。");
            }
        }
        public bool LoadMaterialLibrary(string path)
        {
            _mtllib = new MaterialLibrary(path);

            if (_mtllib != null) return true;
            else return false;
        }
        public bool Load(string path)
        {
            return Load(path, false);
        }
        public bool Load(string path, bool ignoreMtllib)
        {
            FileStream fs;
            StreamReader sr;

            try
            {
                fs = new FileStream(path, FileMode.Open);
                sr = new StreamReader(fs);
            }
            catch (Exception ex)
            {
                throw ex;
            }

            List<Vector3d> plist = new List<Vector3d>();
            List<Vector3d> nlist = new List<Vector3d>();
            List<Vector2d> tclist = new List<Vector2d>();

            List<Facet> flist = new List<Facet>();

            List<string> gnamelist = new List<string>();

            int MaterialID = -1;
            int GroupID = 0;
            string s;
            while ((s = sr.ReadLine()) != null)
            {
                if (s.Length == 0) continue;
                if (s.StartsWith("#")) continue;

                string[] token = s.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                switch (token[0].ToLower())
                {
                    case "v":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("パースに失敗:" + token[1]);
                            if (!double.TryParse(token[2], out y)) throw new Exception("パースに失敗:" + token[2]);
                            if (!double.TryParse(token[3], out z)) throw new Exception("パースに失敗:" + token[3]);
                            plist.Add(new Vector3d(x, y, z));
                            break;
                        }

                    case "vn":
                        {
                            double nx, ny, nz;
                            if (!double.TryParse(token[1], out nx)) throw new Exception("パースに失敗:" + token[1]);
                            if (!double.TryParse(token[2], out ny)) throw new Exception("パースに失敗:" + token[2]);
                            if (!double.TryParse(token[3], out nz)) throw new Exception("パースに失敗:" + token[3]);
                            nlist.Add(new Vector3d(nx, ny, nz));
                            break;
                        }
                    case "vt":
                        {
                            double u, v;
                            if (!double.TryParse(token[1], out u)) throw new Exception("パースに失敗:" + token[1]);
                            if (!double.TryParse(token[2], out v)) throw new Exception("パースに失敗:" + token[2]);
                            tclist.Add(new Vector2d(u, v));
                            break;
                        }
                    case "f":
                        {
                            Facet f = new Facet(token);
                            f.MaterialId = MaterialID;
                            f.GroupId = GroupID;
                            flist.Add(f);
                            break;
                        }
                    case "mtllib":
                        {
                            if (!ignoreMtllib) _mtllib = new MaterialLibrary(_dir + "\\" + token[1]);
                            break;
                        }
                    case "usemtl":
                        {
                            MaterialID = _mtllib.MaterialId(token[1]);
                            break;
                        }
                    case "g":
                        {
                            GroupID = gnamelist.Count;
                            gnamelist.Add(token[1]);
                            break;
                        }
                    default:
                        break;
                }

            }

            // 移し変え
            if (plist.Count != 0) _position = plist.ToArray();
            //{
            //    _position = plist.ToArray();// new Vector3d[plist.Count];
            //    for (int i = 0; i < plist.Count; i++) _position[i] = plist[i];
            //}

            if (nlist.Count != 0) _normal = nlist.ToArray();
            //{
            //    _normal = new Vector3d[nlist.Count];
            //    for (int i = 0; i < nlist.Count; i++) _normal[i] = nlist[i];
            //}

            if (tclist.Count != 0) _texCoord = tclist.ToArray();
            //{
            //    _texCoord = new Vector2d[tclist.Count];
            //    for (int i = 0; i < tclist.Count; i++) _texCoord[i] = tclist[i];
            //}

            if (flist.Count != 0) _facets = flist.ToArray();
            //{
            //    _facets = new Facet[flist.Count];
            //    for (int i = 0; i < flist.Count; i++) _facets[i] = flist[i];
            //}

            if (gnamelist.Count != 0) _groupNames = gnamelist.ToArray();
            else _groupNames = new string[] { "default" };


            sr.Close();
            fs.Close();

            return true;
        }
        public bool LoadBinary(string path)
        {

            FileStream fs;
            BinaryReader br;

            try
            {
                fs = new FileStream(path, FileMode.Open);
                br = new BinaryReader(fs);
            }
            catch (Exception ex)
            {
                throw new Exception(ex.ToString());
            }

            int np = br.ReadInt32();
            int ntc = br.ReadInt32();
            int nnor = br.ReadInt32();
            int nf = br.ReadInt32();

            _position = new Vector3d[np];
            _texCoord = new Vector2d[ntc];
            _normal = new Vector3d[nnor];
            _facets = new Facet[nf];


            for (int i = 0; i < _position.Length; i++)
            {
                double x = br.ReadDouble();
                double y = br.ReadDouble();
                double z = br.ReadDouble();
                _position[i] = new Vector3d(x, y, z);
            }

            for (int i = 0; i < _texCoord.Length; i++)
            {
                double x = br.ReadDouble();
                double y = br.ReadDouble();
                _texCoord[i] = new Vector2d(x, y);
            }

            for (int i = 0; i < _normal.Length; i++)
            {
                double x = br.ReadDouble();
                double y = br.ReadDouble();
                double z = br.ReadDouble();
                _normal[i] = new Vector3d(x, y, z);
            }

            for (int i = 0; i < _facets.Length; i++)
            {
                int ne = br.ReadInt32();
                _facets[i] = new Facet(ne);

                int matid = br.ReadInt32();
                _facets[i].MaterialId = matid;

                for (int j = 0; j < ne; j++)
                {
                    _facets[i].PosIdx[j] = br.ReadInt32();
                    _facets[i].TCIdx[j] = br.ReadInt32();
                    _facets[i].NorIdx[j] = br.ReadInt32();
                }
            }
            br.Close();
            fs.Close();
            return true;
        }
        /// <summary>
        /// Y軸最小値
        /// </summary>
        /// <returns></returns>
        public double Floor()
        {
            double tempMin = double.MaxValue;

            for(int i = 0; i < NumVertex; i++)
            {
                if(Position(i).Y < tempMin) tempMin = Position(i).Y;
            }

            return tempMin;
        }

        public void CenterScale(out Vector3d center, out double scale)
        {
            Vector3d p = new Vector3d(0, 0, 0);

            for (int i = 0; i < NumVertex; i++) p += _position[i];

            p = 1.0 / (double)NumVertex * p;

            center = p;


            double MaxLenSq = double.MinValue;

            for (int i = 0; i < NumVertex; i++)
            {
                double lenSq = Vector3d.LengthSq(_position[i] - center);
                if (MaxLenSq < lenSq) MaxLenSq = lenSq;
            }

            scale = 1.0 / Math.Sqrt(MaxLenSq);
        }

        public void SetupMaterialFacet()
        {
            List<int>[] flist = new List<int>[NumMaterials];
            for (int i = 0; i < flist.Length; i++) flist[i] = new List<int>();

            for (int i = 0; i < NumFacet; i++)
            {
                int mid = _facets[i].MaterialId;
                flist[mid].Add(i);
                // mtllib使ってなかったら落ちるなあ。
            }

            _facetIdPermat = new int[NumMaterials][];
            for (int i = 0; i < NumMaterials; i++)
            {
                _facetIdPermat[i] = new int[flist[i].Count];
                for (int j = 0; j < flist[i].Count; j++)
                {
                    _facetIdPermat[i][j] = flist[i][j];
                }
            }
        }

        public void Triangulate()
        {
            List<Facet> flist = new List<Facet>();
            for (int i = 0; i < NumFacet; i++)
            {
                Facet[] f = _facets[i].Triangulate();

                for (int j = 0; j < f.Length; j++) flist.Add(f[j]);
            }

            // 移し変え
            _facets = new Facet[flist.Count];
            for (int i = 0; i < flist.Count; i++)
            {
                _facets[i] = flist[i];
            }

        }

        /// <summary>
        /// 三角形にのみ有効
        /// </summary>
        public void RemoveZeroSurfaceFacet()
        {
            List<Facet> flist = new List<Facet>();
            for (int i = 0; i < _facets.Length; i++)
            {
                Vector3d OA = Position(i, 1) - Position(i, 0);
                Vector3d OB = Position(i, 2) - Position(i, 0);

                // 面積
                double S = 0.5 * Math.Sqrt(Vector3d.LengthSq(OA) * Vector3d.LengthSq(OB) - Vector3d.Dot(OA, OB) * Vector3d.Dot(OA, OB));

                if (Math.Abs(S) > 0.0001) flist.Add(_facets[i]);
            }

            // 移し替え
            _facets = flist.ToArray();
        }
        Vector3d ComputeTangent(
            Vector3d p0, Vector2d uv0,
            Vector3d p1, Vector2d uv1,
            Vector3d p2, Vector2d uv2)
        {
            // 5次元→3次元頂点に
            Vector3d[] CP0 = new Vector3d[]{
      new Vector3d( p0.X, uv0.X, uv0.Y ),
      new Vector3d( p0.Y, uv0.X, uv0.Y ),
      new Vector3d( p0.Z, uv0.X, uv0.Y ),
   };
            Vector3d[] CP1 = new Vector3d[]{
      new Vector3d( p1.X, uv1.X, uv1.Y ),
      new Vector3d( p1.Y, uv1.X, uv1.Y ),
      new Vector3d( p1.Z, uv1.X, uv1.Y ),
   };
            Vector3d[] CP2 = new Vector3d[]{
      new Vector3d( p2.X, uv2.X, uv2.Y ),
      new Vector3d( p2.Y, uv2.X, uv2.Y ),
      new Vector3d( p2.Z, uv2.X, uv2.Y ),
   };

            // 平面パラメータからUV軸座標算出
            double[] U = new double[3];
            double[] V = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                Vector3d V1 = CP1[i] - CP0[i];
                Vector3d V2 = CP2[i] - CP1[i];
                Vector3d ABC = Vector3d.Cross(V1, V2);

                if (ABC.X == 0.0f)
                {
                    // やばいす！
                    // ポリゴンかUV上のポリゴンが縮退してます！
                    throw new Exception("やばいす！ ポリゴンかUV上のポリゴンが縮退してます！");
                }
                U[i] = -ABC.Y / ABC.X;
                V[i] = -ABC.Z / ABC.X;
            }

            return new Vector3d(U[0], U[1], U[2]);
        }
        public void ComputeTangent()
        {
            // 頂点数と同じだけ接線を設定
            _tangent = new Vector3d[_position.Length];

            // すべての面に対して法線計算し、出てくる頂点にその法線ベクトルを加算
            for (int i = 0; i < NumFacet; i++)
            {
                Vector3d[] vec = new Vector3d[_facets[i].NumEdge];
                Vector2d[] tex = new Vector2d[_facets[i].NumEdge];
                // まず代入
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    vec[j] = Position(i, j);
                    tex[j] = _facets[i].IsTexCoord ? TexCoord(i, j) : new Vector2d();
                }

                Vector3d nor = Vector3d.Cross(Vector3d.Normalize(vec[1] - vec[0]), Vector3d.Normalize(vec[2] - vec[0]));
                Vector3d tan = this.ComputeTangent(vec[0], tex[0], vec[1], tex[1], vec[2], tex[2]);

                double dot = Vector3d.Dot(nor, tan);

                if (Vector3d.IsNaN(tan)) tan = new Vector3d(1, 0, 0);

                // 該当するポリゴンを構成する全頂点に加算
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    _tangent[_facets[i].PosIdx[j]] += Vector3d.Normalize(tan);

                    // 三角形の面積で重みを付けて加算
                    //_tangent[_facets[i].PosIdx[j]] += 0.5 * Math.Sqrt(Vector3d.LengthSq(P) * Vector3d.LengthSq(Q) - Vector3d.Dot(P, Q) * Vector3d.Dot(P, Q)) * Vector3d.Normalize(tan);
                }
            }

            // 接線の正規化
            for (int i = 0; i < _tangent.Length; i++)
            {
                _tangent[i] = Vector3d.Normalize(_tangent[i]);

                if (Vector3d.IsNaN(_tangent[i])) _tangent[i] = new Vector3d(1, 0, 0);
            }


            for (int i = 0; i < NumFacet; i++)
            {
                // TanIdxを作成
                _facets[i].SetupTangentIdx(_facets[i].NumEdge);

                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    // Tanに登録
                    _facets[i].TanIdx[j] = _facets[i].PosIdx[j];
                }
            }
        }
        public void ComputeTangent_()
        {
            // 頂点数と同じだけ接線を設定
            _tangent = new Vector3d[_position.Length];

            // すべての面に対して法線計算し、出てくる頂点にその法線ベクトルを加算
            for (int i = 0; i < NumFacet; i++)
            {
                Vector3d[] vec = new Vector3d[_facets[i].NumEdge];
                Vector2d[] tex = new Vector2d[_facets[i].NumEdge];
                Vector3d tan;

                // まず代入
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    vec[j] = Position(i, j);
                    tex[j] = _facets[i].IsTexCoord ? TexCoord(i, j) : new Vector2d();
                }

                // 接線計算
                Vector3d P = vec[1] - vec[0];
                Vector3d Q = vec[2] - vec[0];

                double s1 = tex[1].X - tex[0].X;
                double t1 = tex[1].Y - tex[0].Y;
                double s2 = tex[2].X - tex[0].X;
                double t2 = tex[2].Y - tex[0].Y;

                double Det = 1.0 / (s1 * t2 - s2 * t1);

                double tx = Det * (t2 * P.X - t1 * Q.X);
                double ty = Det * (t2 * P.Y - t1 * Q.Y);
                double tz = Det * (t2 * P.Z - t1 * Q.Z);

                // 外積計算
                tan = new Vector3d(tx, ty, tz);

                tan = Vector3d.Normalize(tan);

                if (Vector3d.IsNaN(tan)) tan = new Vector3d(1, 0, 0);
   
                // 該当するポリゴンを構成する全頂点に加算
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    //_tangent[_facets[i].PosIdx[j]] += Vector3d.Normalize(tan);

                    // 三角形の面積で重みを付けて加算
                    _tangent[_facets[i].PosIdx[j]] += 0.5 * Math.Sqrt(Vector3d.LengthSq(P) * Vector3d.LengthSq(Q) - Vector3d.Dot(P, Q) * Vector3d.Dot(P, Q)) * Vector3d.Normalize(tan);
                }
            }

            // 接線の正規化
            for (int i = 0; i < _tangent.Length; i++)
            {
                _tangent[i] = Vector3d.Normalize(_tangent[i]);

                if (Vector3d.IsNaN(_tangent[i])) _tangent[i] = new Vector3d(1, 0, 0);
            }


            for (int i = 0; i < NumFacet; i++)
            {
                // TanIdxを作成
                _facets[i].SetupTangentIdx(_facets[i].NumEdge);

                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    // Tanに登録
                    _facets[i].TanIdx[j] = _facets[i].PosIdx[j];
                }
            }

            // binormal = cross(normal, tangent)
        }

        public void SmoothNormal()
        {
            Vector3d[] normals = new Vector3d[_normal.Length];
            for (int i = 0; i < _normal.Length; i++) normals[i] = new Vector3d(0, 0, 0);

            for (int i = 0; i < NumFacet; i++)
            {
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    Vector3d n = Normal(i, j);

                    normals[_facets[i].NorIdx[j]] += n;                    
                }
            }

            for (int i = 0; i < _normal.Length; i++)
            {
                _normal[i] = Vector3d.Normalize(normals[i]);
            }
        }

        public void ComputeNormal()
        {
            // 3角形ポリゴンが前提
            _normal = new Vector3d[_position.Length];

            for (int i = 0; i < NumFacet; i++)
            {
                Vector3d p0 = Position(i, 0);
                Vector3d p1 = Position(i, 1);
                Vector3d p2 = Position(i, 2);

                Vector3d n = Vector3d.Cross(p1 - p0, p2 - p0);
                n.Normalize();

                _normal[_facets[i].PosIdx[0]] = n;
                _normal[_facets[i].PosIdx[1]] = n;
                _normal[_facets[i].PosIdx[2]] = n;

                _facets[i].NorIdx[0] = _facets[i].PosIdx[0];
                _facets[i].NorIdx[1] = _facets[i].PosIdx[1];
                _facets[i].NorIdx[2] = _facets[i].PosIdx[2];
            }

        }

        public void Save(string path)
        {
            FileStream fs;
            StreamWriter sw;

            try
            {
                fs = new FileStream(path, FileMode.Create);
                sw = new StreamWriter(fs, Encoding.GetEncoding("Shift_JIS"));
            }
            catch (Exception ex)
            {
                throw ex;
            }


            for (int i = 0; i < _position.Length; i++)
            {
                Vector3d p = _position[i];
                sw.WriteLine("v " + _position[i].X.ToString() + " " + _position[i].Y.ToString() + " " + _position[i].Z.ToString());
            }


            if (_normal != null)
            {
                for (int i = 0; i < _normal.Length; i++)
                {
                    Vector3d n = _normal[i];

                    sw.WriteLine("vn " + _normal[i].X.ToString() + " " + _normal[i].Y.ToString() + " " + _normal[i].Z.ToString());
                }
            }

            if (_texCoord != null)
            {
                for (int i = 0; i < _texCoord.Length; i++)
                {
                    sw.WriteLine("vt " + _texCoord[i].X.ToString() + " " + _texCoord[i].Y.ToString());
                }
            }

            for (int i = 0; i < _facets.Length; i++)
            {
                sw.Write("f");

                Facet f = _facets[i];
                for (int j = 0; j < f.NumEdge; j++)
                {
                    sw.Write(" ");
                    sw.Write((f.PosIdx[j] + 1).ToString());
                    if (f.IsTexCoord) sw.Write("/" + (f.TCIdx[j] + 1).ToString());
                    if (f.IsNormal) sw.Write("/" + (f.NorIdx[j] + 1).ToString());
                }
                sw.Write("\n");

            }


            sw.Close();
            fs.Close();
        }

        public void Normalize()
        {
            Vector3d center;
            double scale;
            CenterScale(out center, out scale);

            for (int i = 0; i < NumVertex; i++)
            {
                _position[i] = scale * (_position[i] - center);
            }
        }
        public void Scale(double s)
        {
            for (int i = 0; i < NumVertex; i++) _position[i] = s * _position[i];
        }


        public int NumFacet { get { return _facets.Length; } }
        public int NumVertex { get { return _position.Length; } }

        public Vector3d Position(int idx) { return _position[idx]; }
        public Vector3d Position(int fidx, int vidx) { return _position[_facets[fidx].PosIdx[vidx]]; }
        public Vector3d Normal(int idx) { return _normal[idx]; }
        public Vector3d Normal(int fidx, int nidx) { return _normal[_facets[fidx].NorIdx[nidx]]; }
        public Vector2d TexCoord(int idx) { return _texCoord[idx]; }
        public Vector2d TexCoord(int fidx, int tcidx) { return _texCoord[_facets[fidx].TCIdx[tcidx]]; }
        public Vector3d Position(int matid, int fidx, int vidx) { return _position[_facets[_facetIdPermat[matid][fidx]].PosIdx[vidx]]; }
        public Vector3d Normal(int matid, int fidx, int vnidx) { return _normal[_facets[_facetIdPermat[matid][fidx]].NorIdx[vnidx]]; }
        public Vector2d TexCoord(int matid, int fidx, int tcidx) { return _texCoord[_facets[_facetIdPermat[matid][fidx]].TCIdx[tcidx]]; }
        public Vector3d Tangent(int idx) { return _tangent[idx]; }
        public Vector3d Tangent(int fidx, int vidx) { return _tangent[_facets[fidx].TanIdx[vidx]]; }
        public Vector3d Tangent(int matid, int fidx, int vidx) { return _tangent[_facets[_facetIdPermat[matid][fidx]].TanIdx[vidx]]; }

        public Facet GetFacet(int matid, int fidx) { return _facets[_facetIdPermat[matid][fidx]]; } 

        public bool IsTexCoord { get { return _texCoord != null; } }
        public bool IsNormal    { get { return _normal != null; } }
        public bool IsTangent { get { return _tangent != null; } }

        public double MinX
        {
            get
            {
                double min = double.MaxValue;
                foreach (Vector3d p in _position) if (p.X < min) min = p.X;
                return min;
            }
        }

        public double MinY
        {
            get
            {
                double min = double.MaxValue;
                foreach (Vector3d p in _position) if (p.Y < min) min = p.Y;
                return min;
            }
        }
        public double MinZ
        {
            get
            {
                double min = double.MaxValue;
                foreach (Vector3d p in _position) if (p.Z < min) min = p.Z;
                return min;
            }
        }
        public double MaxX
        {
            get
            {
                double max = double.MinValue;
                foreach (Vector3d p in _position) if (p.X > max) max = p.X;
                return max;
            }
        }

        public double MaxY
        {
            get
            {
                double max = double.MinValue;
                foreach (Vector3d p in _position) if (p.Y > max) max = p.Y;
                return max;
            }
        }
        public double MaxZ
        {
            get
            {
                double max = double.MinValue;
                foreach (Vector3d p in _position) if (p.Z > max) max = p.Z;
                return max;
            }
        }

        public Facet[] Facets
        {
            get { return _facets; }
        }

        /// <summary>
        /// インデックスバッファ作成用に。
        /// </summary>
        /// <param name="n"></param>
        public void AlignNormal(out Vector3d[]n)
        {
            n = new Vector3d[_position.Length];

            for (int i = 0; i < NumFacet; i++)
            {
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    int pidx = _facets[i].PosIdx[j];

                    n[pidx] = _normal[_facets[i].NorIdx[j]];
                }
            }

        }

        public void AlignTexCoord(out Vector2d[] tc)
        {
            tc = new Vector2d[_position.Length];

            for (int i = 0; i < NumFacet; i++)
            {
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    int pidx = _facets[i].PosIdx[j];

                    tc[pidx] = _facets[i].IsTexCoord ? _texCoord[_facets[i].TCIdx[j]] : new Vector2d();
                }
            }
        }


        public void SaveAsBin(string path)
        {
            FileStream fs;
            BinaryWriter bw;

            try
            {
                fs = new FileStream(path, FileMode.Create);
                bw = new BinaryWriter(fs);
            }
            catch (Exception ex)
            {
                throw new Exception(ex.ToString());
            }

            bw.Write((Int32)_position.Length);
            bw.Write((Int32)_texCoord.Length);
            bw.Write((Int32)_normal.Length);
            bw.Write((Int32)NumFacet);

            // Position
            for (int i = 0; i < _position.Length; i++)
            {
                bw.Write((double)_position[i].X);
                bw.Write((double)_position[i].Y);
                bw.Write((double)_position[i].Z);
            }

            for (int i = 0; i < _texCoord.Length; i++)
            {
                bw.Write((double)_texCoord[i].X);
                bw.Write((double)_texCoord[i].Y);
            }

            for (int i = 0; i < _normal.Length; i++)
            {
                bw.Write((double)_normal[i].X);
                bw.Write((double)_normal[i].Y);
                bw.Write((double)_normal[i].Z);
            }
            // Facet
            for (int i = 0; i < NumFacet; i++)
            {
                bw.Write((Int32)_facets[i].NumEdge);
                bw.Write((Int32)_facets[i].MaterialId);
                for (int j = 0; j < _facets[i].NumEdge; j++)
                {
                    bw.Write((Int32)_facets[i].PosIdx[j]);
                    bw.Write((Int32)_facets[i].TCIdx[j]);
                    bw.Write((Int32)_facets[i].NorIdx[j]);
                }
            }
                
            



            bw.Close();
            fs.Close();
        }

        public static WavefrontObj FromBin(string path)
        {
            WavefrontObj obj = null;
            using (Stream stream = File.OpenRead(path))
            {
                BinaryFormatter formatter = new BinaryFormatter();

                obj = (WavefrontObj)formatter.Deserialize(stream);
            }

            return obj;
        }
    }

    //[Serializable]
    public class WavefrontMaterial
    {
        public string Name;

        /// <summary>
        /// Shininessだと思う。
        /// </summary>
        public double Ns;

        /// <summary>
        /// 光の屈折率
        /// </summary>
        public double Ni;

        /// <summary>
        /// アルファ
        /// </summary>
        public double d;

        /// <summary>
        /// ????
        /// </summary>
        public double Tr;

        /// <summary>
        /// ????
        /// </summary>
        public Vector3d Tf;

        /// <summary>
        /// 0:照明なし、1:反射ハイライトなし、2:Ksの値で反射ハイライトあり
        /// </summary>
        public int Illum;

        /// <summary>
        /// 環境色
        /// </summary>
        public Vector3d Ka;

        /// <summary>
        /// 拡散反射係数
        /// </summary>
        public Vector3d Kd;

        /// <summary>
        /// スペキュラ係数
        /// </summary>
        public Vector3d Ks;

        /// <summary>
        /// エミッション
        /// </summary>
        public Vector3d Ke;

        public int MapKdId;
        public int MapKaId;
        public int MapBumpId;
        /// <summary>
        /// エフェクト名
        /// </summary>
        public string Technique;

        public double Eta;

        public bool IdxBufferReady;

        public bool Mirror;

        /// <summary>
        /// Shrinkの近似式で使用する。入射角0[deg]のときのフレネル反射率。
        /// </summary>
        public double F0
        {
            get
            {
                return 0.5 * (((Eta - 1.0) * (Eta - 1.0) + (1.0 - Eta) * (1.0 - Eta)) / ((1.0 + Eta) * (1.0 + Eta)));
            }
        }



        public WavefrontMaterial(
            string name, double ns, double ni, double d, double tr,
            Vector3d tf, int illum, Vector3d ka, Vector3d kd, Vector3d ks, Vector3d ke, int mapKdId, int mapBumpId, string tech, double eta, bool idxbuffer, bool mirror)
        {
            this.Name = name;
            this.Ns = ns;
            this.Ni = ni;
            this.d = d;
            this.Tr = tr;
            this.Tf = tf;
            this.Illum = illum;
            this.Ka = ka;
            this.Kd = kd;
            this.Ks = ks;
            this.Ke = ke;
            this.Technique = tech;
            this.MapKdId = mapKdId;
            this.MapBumpId = mapBumpId;
            this.Eta = eta;
            this.IdxBufferReady = idxbuffer;
            this.Mirror = mirror;
        }
        public WavefrontMaterial()
            : this(
            "default", 0.0, 1.5, 1.0, 0.0, new Vector3d(1, 1, 1), 2,
            new Vector3d(0, 0, 0), new Vector3d(1, 1, 1), new Vector3d(0, 0, 0), new Vector3d(0, 0, 0), -1, -1, "default", 1.0, false, false)
        {
        }
    };

    public class MaterialLibrary
    {
        private WavefrontMaterial[] _materials;
        public int NumMaterials { get { return _materials.Length; } }
        public WavefrontMaterial[] Materials { get { return _materials; } }

        private string[] _mapNames; // マップの名前
        public int NumMaps{get{return _mapNames.Length;}}
        public string[] MapNames { get { return _mapNames; } }

        public void SetKdByName(string name, Vector3d kd)
        {
            int id = MaterialId(name);
            if (id != -1)
            {
                _materials[id].Kd = kd;
            }
        }
        public Vector3d GetKdByName(string name)
        {
            int id = MaterialId(name);
            if (id != -1)
            {
                return _materials[id].Kd;
            }
            else throw new Exception("Materialが見つからない。");
        }

        public int MaterialId(string name)
        {
            for (int i = 0; i < _materials.Length; i++)
            {
                if (_materials[i].Name == name) return i;
            }

            throw new Exception("マテリアルが見つからない:" + name);
            return -1;
        }
        public MaterialLibrary()
        {
            _materials = new WavefrontMaterial[1];
            _materials[0] = new WavefrontMaterial();
        }
        public MaterialLibrary(string path)
        {
            if (!Load(path))
            {
                return;
            }
        }

        public bool Load(string path)
        {
            FileStream fs;
            StreamReader sr;

            try
            {
                fs = new FileStream(path, FileMode.Open);
                sr = new StreamReader(fs);
            }
            catch (Exception ex)
            {
                throw ex;
            }


            WavefrontMaterial mtl = new WavefrontMaterial();
            List<WavefrontMaterial> list = new List<WavefrontMaterial>();
            List<string> maplist = new List<string>();

            int cnt = 0;
            string s;
            while ((s = sr.ReadLine()) != null)
            {
                if (s.Length == 0) continue;
                if (s.StartsWith("#")) continue;

                string[] token = s.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                switch (token[0].ToLower())
                {
                    case "newmtl":
                        {
                            if (cnt != 0) list.Add(mtl);
                            mtl = new WavefrontMaterial();
                            mtl.Name = token[1];
                            cnt++;
                        }
                        break;
                    case "ns":
                        {
                            double ns;
                            if (!double.TryParse(token[1], out ns)) throw new Exception("不正なフォーマット");
                            else mtl.Ns = ns;
                        }
                        break;
                    case "ni":
                        {
                            double ni;
                            if (!double.TryParse(token[1], out ni)) throw new Exception("不正なフォーマット");
                            else mtl.Ni = ni;
                        }
                        break;
                    case "d":
                        {
                            double d;
                            if (!double.TryParse(token[1], out d)) throw new Exception("不正なフォーマット");
                            else mtl.d = d;
                        }
                        break;
                    case "tr":
                        {
                            double tr;
                            if (!double.TryParse(token[1], out tr)) throw new Exception("不正なフォーマット");
                            else mtl.Tr = tr;
                        }
                        break;
                    case "tf":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[2], out y)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[3], out z)) throw new Exception("不正なフォーマット");
                            else mtl.Tf = new Vector3d(x, y, z);
                        }
                        break;
                    case "illum":
                        {
                            int illum;
                            if (!int.TryParse(token[1], out illum)) throw new Exception("不正なフォーマット");
                            else mtl.Illum = illum;
                        }
                        break;
                    case "ka":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[2], out y)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[3], out z)) throw new Exception("不正なフォーマット");
                            else mtl.Ka = new Vector3d(x, y, z);
                        }
                        break;
                    case "kd":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[2], out y)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[3], out z)) throw new Exception("不正なフォーマット");
                            else mtl.Kd = new Vector3d(x, y, z);
                        }
                        break;
                    case "ks":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[2], out y)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[3], out z)) throw new Exception("不正なフォーマット");
                            else mtl.Ks = new Vector3d(x, y, z);
                        }
                        break;
                    case "ke":
                        {
                            double x, y, z;
                            if (!double.TryParse(token[1], out x)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[2], out y)) throw new Exception("不正なフォーマット");
                            if (!double.TryParse(token[3], out z)) throw new Exception("不正なフォーマット");
                            else mtl.Ke = new Vector3d(x, y, z);
                        }
                        break;
                    case "map_kd":
                        {
                            int id;
                            if (-1 == (id = maplist.IndexOf(token[1])))
                            {
                                mtl.MapKdId = maplist.Count;
                                maplist.Add(token[1]);
                            }
                            else
                            {
                                mtl.MapKdId = id;
                            }
                            break;
                        }
                    case "map_ka":
                        {
                            int id;
                            if (-1 == (id = maplist.IndexOf(token[1])))
                            {
                                mtl.MapKaId = maplist.Count;
                                maplist.Add(token[1]);
                            }
                            else
                            {
                                mtl.MapKaId = id;
                            }
                            break;
                        }
                    case "map_bump":
                        {
                            int id;
                            if (-1 == (id = maplist.IndexOf(token[1])))
                            {
                                mtl.MapBumpId = maplist.Count;
                                maplist.Add(token[1]);
                            }
                            else
                            {
                                mtl.MapBumpId = id;
                            }
                            break;
                        }
                    case "technique":
                        {
                            mtl.Technique = token[1];
                        }
                        break;
                    case "eta": // 空気との相対屈折率
                        {
                            mtl.Eta = double.Parse(token[1]);
                        }
                        break;
                    case "idxbuffer":
                        {
                            mtl.IdxBufferReady = bool.Parse(token[1]);
                            break;
                        }
                    case "mirror":
                        {
                            mtl.Mirror = bool.Parse(token[1]);
                            break;
                        }
                    default:
                        System.Diagnostics.Debug.WriteLine("不明なトークン" + token[0]);
                        break;
                }

            }

            // 最後にも追加
            list.Add(mtl);

            // 移し変え
            _materials = new WavefrontMaterial[list.Count];
            for (int i = 0; i < list.Count; i++) _materials[i] = list[i];

            _mapNames = new string[maplist.Count];
            for (int i = 0; i < maplist.Count; i++) _mapNames[i] = maplist[i];

            sr.Close();
            fs.Close();

            return true;
        }
        public void TransparentMaterialSort()
        {
            //並び替える配列を作成
            double[] trans = new double[_materials.Length];
            WavefrontMaterial[] mtls = new WavefrontMaterial[_materials.Length];
            for (int i = 0; i < _materials.Length; i++)
            {
                trans[i] = _materials[i].d;
                mtls[i] = _materials[i];
            }

            Array.Sort(trans, mtls);

            // 降順
            for (int i = 0; i < mtls.Length; i++)
            {
                _materials[i] = mtls[mtls.Length - i - 1];
            }

        }
    };
}
