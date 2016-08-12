using System;
using System.Collections.Generic;
using System.Text;

namespace MathUtil
{
    [System.Diagnostics.DebuggerDisplay("X:{_x}, Y:{_y}, Z:{_z}, W:{_w}")]
    public struct Vector4d
    {
        private double _x, _y, _z, _w;

        public double X
        {
            get { return _x; }
            set { _x = value; }
        }
        public double Y
        {
            get { return _y; }
            set { _y = value; }
        }
        public double Z
        {
            get { return _z; }
            set { _z = value; }
        }
        public double W
        {
            get { return _w; }
            set { _w = value; }
        }


        /// <summary>
        /// 配列風にアクセスするインデクサ
        /// </summary>
        /// <param name="idx"></param>
        /// <returns></returns>
        public double this[int idx]
        {
            get
            {
                if (idx == 0) return _x;
                else if (idx == 1) return _y;
                else if (idx == 2) return _z;
                else if( idx == 3)  return _w;
                else throw new Exception("indexが範囲外です。");
            }
            set
            {
                if (idx == 0) _x = value;
                else if (idx == 1) _y = value;
                else if (idx == 2) _z = value;
                else if (idx == 3) _w = value;
                else throw new Exception("indexが範囲外です。");
            }
        }

        public Vector4d(double x, double y, double z, double w)
        {
            _x = x;
            _y = y;
            _z = z;
            _w = w;
        }
        public Vector4d(Vector3d v, double w)
        {
            _x = v.X;
            _y = v.Y;
            _z = v.Z;
            _w = w;
        }

        public static Vector4d FromSphereCoord(double len, double theta, double phi)
        {
            double x = Math.Sin(theta) * Math.Cos(phi);
            double y = Math.Sin(theta) * Math.Sin(phi);
            double z = Math.Cos(theta);

            return len * new Vector4d(x, y, z, 1);
        }
        public static double Dot(Vector4d v1, Vector4d v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z + v1.W * v2.W;
        }
        public static double LengthSq(Vector4d v)
        {
            return v.X * v.X + v.Y * v.Y + v.Z * v.Z + v.W * v.W;
        }
        public static double Length(Vector4d v)
        {
            double LenSq = Vector4d.LengthSq(v);
            return Math.Sqrt(LenSq);
        }
        public static Vector4d Normalize(Vector4d v)
        {
            double len = Length(v);
            return 1.0 / len * v;
        }
        public static Vector4d operator *(double left, Vector4d right)
        {
            return new Vector4d(left * right.X, left * right.Y, left * right.Z, left * right.W);
        }
        public static Vector4d operator *(Vector4d left, double right)
        {
            return new Vector4d(right * left.X, right * left.Y, right * left.Z, right * left.W);
        }
        public static Vector4d operator +(Vector4d left, Vector4d right)
        {
            return new Vector4d(left.X + right.X, left.Y + right.Y, left.Z + right.Z, left.W + right.W);
        }
        public static Vector4d operator +(Vector4d left, double right)
        {
            return new Vector4d(left.X + right, left.Y + right, left.Z + right, left.W + right);
        }
        public static Vector4d operator +(double left, Vector4d right)
        {
            return new Vector4d(left + right.X, left + right.Y, left + right.Z, left + right.W);
        }
        public static Vector4d operator -(Vector4d left, double right)
        {
            return new Vector4d(left.X - right, left.Y - right, left.Z - right, left.W - right);
        }
        public static Vector4d operator -(double left, Vector4d right)
        {
            return new Vector4d(left - right.X, left - right.Y, left - right.Z, left - right.W);
        }
        public static Vector4d operator -(Vector4d left, Vector4d right)
        {
            return new Vector4d(left.X - right.X, left.Y - right.Y, left.Z - right.Z, left.W - right.W);
        }
        public static Vector4d operator *(Vector4d left, Vector4d right)
        {
            return new Vector4d(left.X * right.X, left.Y * right.Y, left.Z * right.Z, left.W * right.W);
        }
        public static Vector4d operator /(Vector4d left, Vector4d right)
        {
            return new Vector4d(left.X / right.X, left.Y / right.Y, left.Z / right.Z, left.W / right.W);
        }
        // 単項演算子
        public static Vector4d operator -(Vector4d v)
        {
            return new Vector4d(-v.X, -v.Y, -v.Z, -v.W);
        }

        public void Normalize()
        {
            double len = Vector4d.Length(this);
            this = 1.0 / len * this;
        }
   
       
        public static bool IsNaN(Vector4d v)
        {
            return double.IsNaN(v.X) || double.IsNaN(v.Y) || double.IsNaN(v.Z) || double.IsNaN(v.W);
        }

        public Vector3d ToVector3d()
        {
            return new Vector3d(this.X, this.Y, this.Z);
        }
        public void Write(System.IO.BinaryWriter bw)
        {
            bw.Write(this.X);
            bw.Write(this.Y);
            bw.Write(this.Z);
            bw.Write(this.W);
        }
        public static Vector4d FromBinaryReader(System.IO.BinaryReader br)
        {
            return new Vector4d(
            br.ReadDouble(),
            br.ReadDouble(),
            br.ReadDouble(),
            br.ReadDouble());
        }
    }
}
