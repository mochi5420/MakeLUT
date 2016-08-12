using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class RandomVector
    {
        // uに正の方向の半球状に向かうランダムなベクトル
        static public Vector3f RandomUnitVectorHemisphere(Vector3f u, Random rnd)
        {
            Vector3f rv = RandomUnitVector(rnd);

            if (Vector3f.Dot(u, rv) < 0) return -rv;
            else return rv;
        }

        /// <summary>
        /// ランダムな方向の３次元単位ベクトル
        /// </summary>
        /// <param name="rnd"></param>
        /// <returns></returns> 
        static public Vector3f RandomUnitVector(Random rnd)
        {
            Vector3f v;

            while (true)
            {
                // [-1, +1]の乱数を各成分にもつ
                v = new Vector3f((float)(2.0 * rnd.NextDouble() - 1.0), (float)(2.0 * rnd.NextDouble() - 1.0), (float)(2.0 * rnd.NextDouble() - 1.0));


                if (v.LengthSq() > 1.0) continue;
                else break;
            }
            v.Normalize();

            return v;
        }
    }
}
