using System;
using System.Collections.Generic;
using System.Text;

namespace MathUtil
{
    //[Serializable]
    [System.Diagnostics.DebuggerDisplay("X:{_x}, Y:{_y}, Z:{_z}")]
    public struct Vector3d
    {
        private double _x, _y, _z;

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

        public double R
        {
            get { return _x; }
            set { _x = value; }
        }
        public double G
        {
            get { return _y; }
            set { _y = value; }
        }
        public double B
        {
            get { return _z; }
            set { _z = value; }
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
                else throw new Exception("indexが範囲外です。");
            }
            set
            {
                if (idx == 0) _x = value;
                else if (idx == 1) _y = value;
                else if (idx == 2) _z = value;
                else throw new Exception("indexが範囲外です。");
            }
        }

        public Vector3d(double x, double y, double z)
        {
            _x = x;
            _y = y;
            _z = z;
        }
        public override string ToString()
        {
            return "(" + _x.ToString() + ", " + _y.ToString() + ", " + _z.ToString() + ")";
        }
        public static Vector3d FromSphereCoord(double len, double theta, double phi)
        {
            double x = Math.Sin(theta) * Math.Cos(phi);
            double y = Math.Sin(theta) * Math.Sin(phi);
            double z = Math.Cos(theta);

            return len * new Vector3d(x, y, z);
        }
        public static Vector3d FromSphereCoord(double len, double theta, double phi, Vector3d e1, Vector3d e2, Vector3d e3)
        {
            return len * (Math.Sin(theta) * Math.Cos(phi) * e1 + Math.Sin(theta) * Math.Sin(phi) * e2 + Math.Cos(theta) * e3);
        }
        public static double Dot(Vector3d v1, Vector3d v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }
        public static double LengthSq(Vector3d v)
        {
            return v.X * v.X + v.Y * v.Y + v.Z * v.Z;
        }
        public static double Length(Vector3d v)
        {
            double LenSq = Vector3d.LengthSq(v);
            return Math.Sqrt(LenSq);
        }
        public static Vector3d Normalize(Vector3d v)
        {
            double len = Length(v);
            return 1.0 / len * v;
        }
        public static Vector3d operator *(double left, Vector3d right)
        {
            return new Vector3d(left * right.X, left * right.Y, left * right.Z);
        }
        public static Vector3d operator *(Vector3d left, double right)
        {
            return new Vector3d(right * left.X, right * left.Y, right * left.Z);
        }
        public static Vector3d operator +(Vector3d left, Vector3d right)
        {
            return new Vector3d(left.X + right.X, left.Y + right.Y, left.Z + right.Z);
        }
        public static Vector3d operator +(Vector3d left, double right)
        {
            return new Vector3d(left.X + right, left.Y + right, left.Z + right);
        }
        public static Vector3d operator +(double left, Vector3d right)
        {
            return new Vector3d(left + right.X, left + right.Y, left + right.Z);
        }
        public static Vector3d operator -(Vector3d left, double right)
        {
            return new Vector3d(left.X - right, left.Y - right, left.Z - right);
        }
        public static Vector3d operator -(double left, Vector3d right)
        {
            return new Vector3d(left - right.X, left - right.Y, left - right.Z);
        }
        public static Vector3d operator -(Vector3d left, Vector3d right)
        {
            return new Vector3d(left.X - right.X, left.Y - right.Y, left.Z - right.Z);
        }
        public static Vector3d operator *(Vector3d left, Vector3d right)
        {
            return new Vector3d(left.X * right.X, left.Y * right.Y, left.Z * right.Z);
        }
        public static Vector3d operator /(Vector3d left, Vector3d right)
        {
            return new Vector3d(left.X / right.X, left.Y / right.Y, left.Z / right.Z);
        }
        public static Vector3d operator /(Vector3d left, double right)
        {
            return new Vector3d(left.X / right, left.Y / right, left.Z / right);
        }
        // 単項演算子
        public static Vector3d operator -(Vector3d v)
        {
            return new Vector3d(-v.X, -v.Y, -v.Z);
        }
        public static Vector3d Cross(Vector3d v1, Vector3d v2)
        {
            double x = v1.Y * v2.Z - v1.Z * v2.Y;
            double y = v1.Z * v2.X - v1.X * v2.Z;
            double z = v1.X * v2.Y - v1.Y * v2.X;

            return new Vector3d(x, y, z);
        }
        public void Normalize()
        {
            double len = Vector3d.Length(this);
            this = 1.0 / len * this;
        }
        public static Vector3d RandomUnitVector(Random rnd)
        {           

            while (true)
            {
                Vector3d p = new Vector3d(2.0 * rnd.NextDouble() - 1.0, 2.0 * rnd.NextDouble() - 1.0, 2.0 * rnd.NextDouble() - 1.0);

                if (LengthSq(p) < 1.0)
                {
                    p.Normalize();

                    return p;
                }               
            }
        }
        public static Vector3d RandomUnitVectorOnHemiSphere(Vector3d up, Random rnd)
        {
            while (true)
            {
                Vector3d p = RandomUnitVector(rnd);

                if (Vector3d.Dot(p, up) >= 0.0)
                {
                    return p;
                }
            }
        }
        public static Vector3d Rotate(Vector3d source, double theta, double phi)
        {
            Vector3d e1 = source;
            Vector3d e2 = new Vector3d(0, 0, 0);
            Vector3d e3 = new Vector3d(0, 0, 0);

            if (e1.Y == 0.0 && e1.Z == 0.0)
            {
                e2.Y = (double)Math.Sign(e1.X);
                e2.X = e2.Z = 0.0;
            }
            else
            {
                e2.X = 0.0;
                e2.Y = e1.Z;
                e2.Z = -e1.Y;
                e2.Normalize();
            }

            // make e3
            e3 = Vector3d.Cross(e1, e2);

            Vector3d e = Math.Cos(phi) * e2 + Math.Sin(phi) * e3;

            Vector3d vec = Math.Cos(theta) * source + Math.Sin(theta) * e;

            return vec;
        }
        public void ToSphereCoord(out double len, out double theta, out double phi)
        {
            len = Vector3d.Length(this);

            double x = _x / len;
            double y = _y / len;
            double z = _z / len;

            theta = Math.Atan2(Math.Sqrt(x * x + y * y), z);
            phi = Math.Atan2(y, x);
        }

        public static bool IsNaN(Vector3d v)
        {
            return double.IsNaN(v.X) || double.IsNaN(v.Y) || double.IsNaN(v.Z); 
        }
        public static Vector3d Parse(string s, char[]separator)
        {
            string[] token = s.Split(separator);
            double x = double.Parse(token[0]);
            double y = double.Parse(token[1]);
            double z = double.Parse(token[2]);

            return new Vector3d(x, y, z);
        }
        public void Write(System.IO.BinaryWriter bw)
        {
            bw.Write(this.X);
            bw.Write(this.Y);
            bw.Write(this.Z);
        }
        public static Vector3d FromBinaryReader(System.IO.BinaryReader br)
        {
            return new Vector3d(
            br.ReadDouble(),
            br.ReadDouble(),
            br.ReadDouble());
        }
    }
}
