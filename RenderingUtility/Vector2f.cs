using System;

namespace Utility
{
    [System.Diagnostics.DebuggerDisplay("X{X}, Y{Y}")]
    public struct Vector2f
    {
        public float X;
        public float Y;

        public Vector2f(float x, float y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// アクセスが遅いのであまり多用しないこと。
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        public float this[int i]
        {
            set
            {
                if (i == 0) this.X = value;
                else if (i == 1) this.Y = value;
                else
                {
                    throw new Exception("Unsupported Index.");
                }
            }
            get
            {
                if (i == 0) return this.X;
                else if (i == 1) return this.Y;
                else
                {
                    throw new Exception("Unsupported Index.");
                }
            }
        }
        public float Length()
        {
            return Vector2f.Length(this);
        }
        public float LengthSq()
        {
            return Vector2f.LengthSq(this);
        }
        public void Normalize()
        {
            float len = this.Length();
            X = X / len;
            Y = Y / len;
        }

        public static Vector2f operator +(Vector2f left, Vector2f right)
        {
            return new Vector2f(left.X + right.X, left.Y + right.Y);
        }
        public static Vector2f operator -(Vector2f left, Vector2f right)
        {
            return new Vector2f(left.X - right.X, left.Y - right.Y);
        }
        public static Vector2f operator -(Vector2f source)
        {
            return new Vector2f(-source.X, -source.Y);
        }
        public static Vector2f operator *(float left, Vector2f right)
        {
            return new Vector2f(left * right.X, left * right.Y);
        }
        public static float Dot(Vector2f left, Vector2f right)
        {
            return left.X * right.X + left.Y * right.Y;
        }

        public static float Length(Vector2f source)
        {
            return (float)Math.Sqrt(Dot(source, source));
        }
        public static float LengthSq(Vector2f source)
        {
            return Dot(source, source);
        }
        public static Vector2f Normalize(Vector2f source)
        {
            float len = Length(source);
            return 1.0f / len * source;
        }
        public override string ToString()
        {
            string res = "(X:" + X.ToString() + ", Y:" + Y.ToString() + ")";
            return res;
        }
    }
}
