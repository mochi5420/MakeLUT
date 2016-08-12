using System;
using System.Collections.Generic;
using System.Text;

namespace Utility
{
    public struct Matrix33d
    {
        public float M11, M12, M13;
        public float M21, M22, M23;
        public float M31, M32, M33;

        public Matrix33d(float val)
        {
            M11 = M12 = M13 = val;
            M21 = M22 = M23 = val;
            M31 = M32 = M33 = val;
        }

        public Matrix33d(Matrix33d m)
        {
            this = m;
        }


        /// <summary>
        /// アクセスが遅いのであまり多用しないこと。
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public float this[int i, int j]
        {
            get
            {
                if (i == 0 && j == 0) return M11;
                else if (i == 0 && j == 1) return M12;
                else if (i == 0 && j == 2) return M13;
                else if (i == 1 && j == 0) return M21;
                else if (i == 1 && j == 1) return M22;
                else if (i == 1 && j == 2) return M23;
                else if (i == 2 && j == 0) return M31;
                else if (i == 2 && j == 1) return M32;
                else if (i == 2 && j == 2) return M33;
                else throw new Exception("範囲外です。");
            }
            set
            {
                if (i == 0 && j == 0) M11 = value;
                else if (i == 0 && j == 1) M12 = value;
                else if (i == 0 && j == 2) M13 = value;
                else if (i == 1 && j == 0) M21 = value;
                else if (i == 1 && j == 1) M22 = value;
                else if (i == 1 && j == 2) M23 = value;
                else if (i == 2 && j == 0) M31 = value;
                else if (i == 2 && j == 1) M32 = value;
                else if (i == 2 && j == 2) M33 = value;
            }
        }
 
        public static Matrix33d Identity
        {
            get
            {
                Matrix33d I;
                I.M11 = I.M22 = I.M33 = 1.0f;
                I.M12 = I.M13 = I.M21 = I.M23 = I.M31 = I.M32 = 0.0f;
                return I;
            }
        }
        public static Matrix33d Zero
        {
            get
            {
                Matrix33d Zero;
                Zero.M11 = Zero.M12 = Zero.M13 = 0.0f;
                Zero.M21 = Zero.M22 = Zero.M23 = 0.0f;
                Zero.M31 = Zero.M32 = Zero.M33 = 0.0f;
                return Zero;
            }
        }
      
        /// <summary>
        /// 任意軸回転
        /// </summary>
        /// <param name="axisRotation">回転中心になる単位ベクトル。</param>
        /// <param name="angle">回転角（ラジアン）</param>
        /// <returns></returns>
        public static Matrix33d RotateAxis(Vector3f axisRotation, float angle)
        {
            Matrix33d rot = new Matrix33d();
            float Vx,Vy,Vz;
            float cos,sin;
            
            Vx = axisRotation.X;
            Vy = axisRotation.Y;
            Vz = axisRotation.Z;

            cos = (float)Math.Cos(angle);
            sin = (float)Math.Sin(angle);

            rot.M11 = Vx * Vx * (1.0f - cos) + cos; rot.M12 = Vx * Vy * (1.0f - cos) - Vz * sin; rot.M13 = Vz * Vx * (1.0f - cos) + Vy * sin;
            rot.M21 = Vx * Vy * (1.0f - cos) + Vz * sin; rot.M22 = Vy * Vy * (1.0f - cos) + cos; rot.M23 = Vy * Vz * (1.0f - cos) - Vx * sin;
            rot.M31 = Vz * Vx * (1.0f - cos) - Vy * sin; rot.M32 = Vy * Vz * (1.0f - cos) + Vx * sin; rot.M33 = Vz * Vz * (1.0f - cos) + cos;

            return rot;
        }

    }
}
