using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public struct VertexAttribute
    {
        private Vector3f _n;
        private Vector2f _tc;
        private Vector3f _tan;

        private object _object;

        public Vector3f Normal
        {
            get { return _n; }
            set { _n = value; }
        }
        public Vector2f TexCoord
        {
            get { return _tc; }
            set { _tc = value; }
        }
        public Vector3f Tangent
        {
            get { return _tan; }
            set { _tan = value; }
        }
        public object UserObject
        {
            get { return _object; }
            set { _object = value; }
        }


        public VertexAttribute(Vector2f tc, Vector3f n, Vector3f tan)
            : this(tc, n, tan, null)
        {
        }
        public VertexAttribute(Vector2f tc, Vector3f n, Vector3f tan, object userobject)
        {
            _tc = tc;
            _n = n;
            _tan = tan;
            _object = userobject;
        }
    };

    public class Triangle
    {
        private Vector3f[] _v = new Vector3f[3];
       
        private Vector3f _faceNormal;   // ñ ÇÃñ@ê¸


        private VertexAttribute[] _va = new VertexAttribute[3];


        public Triangle(Vector3f[] v)
            : this(v[0], v[1], v[2])
        {

        }
        public Triangle(Vector3f v0, Vector3f v1, Vector3f v2)
            : this(v0, v1, v2, new VertexAttribute(), new VertexAttribute(), new VertexAttribute())
        {
        }
        public Triangle(Vector3f v0, Vector3f v1, Vector3f v2, Vector3f n0, Vector3f n1, Vector3f n2)
            : this(v0, v1, v2,
            new VertexAttribute(new Vector2f(0, 0), n0, new Vector3f(0, 0, 0)),
            new VertexAttribute(new Vector2f(0, 0), n1, new Vector3f(0, 0, 0)),
            new VertexAttribute(new Vector2f(0, 0), n2, new Vector3f(0, 0, 0)))
        {
        }
        public Triangle(Vector3f v0, Vector3f v1, Vector3f v2, VertexAttribute va0, VertexAttribute va1, VertexAttribute va2)
        {
            _v[0] = v0;
            _v[1] = v1;
            _v[2] = v2;

            _va[0] = va0;
            _va[1] = va1;
            _va[2] = va2;

            _faceNormal = Vector3f.Normalize(Vector3f.Cross(_v[1] - _v[0], _v[2] - _v[0]));
        }
       
        public Vector3f V0
        {
            get { return _v[0]; }
            set { _v[0] = value; }
        }
        public Vector3f V1
        {
            get { return _v[1]; }
            set { _v[1] = value; }
        }
        public Vector3f V2
        {
            get { return _v[2]; }
            set { _v[2] = value; }
        }
        public Vector3f[] V
        {
            get { return _v; }
            set { _v = value; }
        }
        public VertexAttribute VA0
        {
            get { return _va[0]; }
            set { _va[0] = value; }
        }
        public VertexAttribute VA1
        {
            get { return _va[1]; }
            set { _va[1] = value; }
        }
        public VertexAttribute VA2
        {
            get { return _va[2]; }
            set { _va[2] = value; }
        }
        public VertexAttribute[] VA
        {
            get { return _va; }
            set { _va = value; }
        }
        public Vector3f FaceNormal
        {
            get { return _faceNormal; }
        }

        public void Thicking()
        {
            this._v[0] += new Vector3f(float.Epsilon, float.Epsilon, float.Epsilon);
            this._v[1] += new Vector3f(float.Epsilon, float.Epsilon, float.Epsilon);
            this._v[2] += new Vector3f(float.Epsilon, float.Epsilon, float.Epsilon);
        }
    }
}
