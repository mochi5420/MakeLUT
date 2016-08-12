using System;

namespace Utility
{
    [System.Diagnostics.DebuggerDisplay("X{X}, Y{Y}, Z{Z}")]
    public struct Vector3f
    {
        public float X;
        public float Y;
        public float Z;

        public Vector3f(float x, float y, float z)
        {
            X = x;
            Y = y;
            Z = z;
        }
        public float Length()
        {
            return Vector3f.Length(this);
        }
        public float LengthSq()
        {
            return Vector3f.LengthSq(this);
        }
        public void Normalize()
        {
            float len = this.Length();
            X = X / len;
            Y = Y / len;
            Z = Z / len;
        }
        public float this[int i]
        {
            set
            {
                if (i == 0) this.X = value;
                else if (i == 1) this.Y = value;
                else if (i == 2) this.Z = value;
                else
                {
                    throw new Exception("Unsupported Index.");
                }
            }
            get
            {
                if (i == 0) return this.X;
                else if (i == 1) return this.Y;
                else if (i == 2) return this.Z;
                else
                {
                    throw new Exception("Unsupported Index.");
                }
            }
        }

        public static Vector3f operator +(Vector3f left, Vector3f right)
        {
            return new Vector3f(left.X + right.X, left.Y + right.Y, left.Z + right.Z);
        }
        public static Vector3f operator -(Vector3f left, Vector3f right)
        {
            return new Vector3f(left.X - right.X, left.Y - right.Y, left.Z - right.Z);
        }
        public static Vector3f operator -(Vector3f source)
        {
            return new Vector3f(-source.X, -source.Y, -source.Z);
        }
        public static Vector3f operator *(float left, Vector3f right)
        {
            return new Vector3f(left * right.X, left * right.Y, left * right.Z);
        }
        public static Vector3f operator *(Vector3f left, Vector3f right)
        {
            return new Vector3f(left.X * right.X, left.Y * right.Y, left.Z * right.Z);
        }
        public static Vector3f operator /(Vector3f left, Vector3f right)
        {
            return new Vector3f(left.X / right.X, left.Y / right.Y, left.Z / right.Z);
        }
        public static float Dot(Vector3f left, Vector3f right)
        {
            return left.X * right.X + left.Y * right.Y + left.Z * right.Z;
        }

        public static Vector3f Cross(Vector3f left, Vector3f right)
        {
            float x = left.Y * right.Z - left.Z * right.Y;
            float y = left.Z * right.X - left.X * right.Z;
            float z = left.X * right.Y - left.Y * right.X;

            return new Vector3f(x, y, z);
        }

        public static float Length(Vector3f source)
        {
            return (float)Math.Sqrt(Dot(source, source));
        }
        public static float LengthSq(Vector3f source)
        {
            return Dot(source, source);
        }
        public static Vector3f Normalize(Vector3f source)
        {
            float len = Length(source);
            return 1.0f / len * source;
        }
        public static Vector3f Mult(Matrix33d sourceMatrix, Vector3f source)
        {
            Vector3f v = new Vector3f(0, 0, 0);
            v.X = sourceMatrix.M11 * source.X + sourceMatrix.M12 * source.Y + sourceMatrix.M13 * source.Z;
            v.Y = sourceMatrix.M21 * source.X + sourceMatrix.M22 * source.Y + sourceMatrix.M23 * source.Z;
            v.Z = sourceMatrix.M31 * source.X + sourceMatrix.M32 * source.Y + sourceMatrix.M33 * source.Z;

            return v;
        }
        public static Vector3f Transform(Vector3f source, Matrix33d sourceMatrix)
        {
            return Mult(sourceMatrix, source);
        }
        public override string ToString()
        {
            string res = "(X:" + X.ToString() + ", Y:" + Y.ToString() + ", Z:" + Z.ToString() + ")";
            return res;
        }
    }
}
