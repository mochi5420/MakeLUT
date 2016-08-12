using System;
using System.Collections.Generic;
using System.Text;

namespace MathUtil
{

//    [Serializable]
    [System.Diagnostics.DebuggerDisplay("X:{_x}, Y:{_y}")]
    public struct Vector2d
    {
        private double _x, _y;

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
                else throw new Exception("indexが範囲外です。");
            }
            set
            {
                if (idx == 0) _x = value;
                else if (idx == 1) _y = value;
                else throw new Exception("indexが範囲外です。");
            }
        }

        public Vector2d(double x, double y)
        {
            _x = x;
            _y = y;
        }
        public override string ToString()
        {
            return "(" + _x.ToString() + ", " + _y.ToString() + ")";
        }
        public static double Dot(Vector2d v1, Vector2d v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y;
        }
        public static double LengthSq(Vector2d v)
        {
            return v.X * v.X + v.Y * v.Y;
        }
        public static double Length(Vector2d v)
        {
            double LenSq = Vector2d.LengthSq(v);
            return Math.Sqrt(LenSq);
        }
        public static Vector2d Normalize(Vector2d v)
        {
            double len = Length(v);
            return 1.0 / len * v;
        }
        public static Vector2d operator *(double left, Vector2d right)
        {
            return new Vector2d(left * right.X, left * right.Y);
        }
        public static Vector2d operator *(Vector2d left, double right)
        {
            return new Vector2d(right * left.X, right * left.Y);
        }
        public static Vector2d operator +(Vector2d left, Vector2d right)
        {
            return new Vector2d(left.X + right.X, left.Y + right.Y);
        }
        public static Vector2d operator -(Vector2d left, Vector2d right)
        {
            return new Vector2d(left.X - right.X, left.Y - right.Y);
        }
        // 単項演算子
        public static Vector2d operator -(Vector2d v)
        {
            return new Vector2d(-v.X, -v.Y);
        }
        public void Noralize()
        {
            double len = Vector2d.Length(this);
            this = 1.0 / len * this;
        }
      

        public void Write(System.IO.BinaryWriter bw)
        {
            bw.Write(this.X);
            bw.Write(this.Y);
        }
        public static Vector2d FromBinaryReader(System.IO.BinaryReader br)
        {
            return new Vector2d(
            br.ReadDouble(),
            br.ReadDouble());
        }
    }
}
