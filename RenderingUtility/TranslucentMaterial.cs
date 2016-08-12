using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class TranslucentMaterial
    {
        /// <summary>
        /// Relative index of refraction
        /// </summary>
        private float _eta;
        /// <summary>
        /// Relative index of refraction
        /// </summary>
        public float Eta
        {
            get { return _eta; }
            set { _eta = value; }
        }

        /// <summary>
        /// Phase function
        /// </summary>
        private float _g;
        /// <summary>
        /// Phase function
        /// </summary>
        public float G
        {
            get { return _g; }
            set { _g = value; }
        }

        /// <summary>
        /// Scattering coefficient
        /// </summary>
        private Vector3d _sigma_s;
        /// <summary>
        /// Scattering coefficient
        /// </summary>
        public Vector3d Sigma_s
        {
            get { return _sigma_s; }
            set { _sigma_s = value; }
        }
        /// <summary>
        /// Extinction coefficient
        /// </summary>
        private Vector3d _sigma_t;
        /// <summary>
        /// Extinction coefficient
        /// </summary>
        public Vector3d Sigma_t
        {
            get { return _sigma_t; }
            set { _sigma_t = value; }
        }

        /// <summary>
        /// scattering albedo
        /// </summary>
        private Vector3d _rambda;
        /// <summary>
        /// scattering albedo
        /// </summary>
        public Vector3d Rambda
        {
            get { return _rambda; }
            set { _rambda = value; }
        }

        /// <summary>
        /// // 垂直入射におけるフレネル反射係数の実部[フォトンマップ本 pp. 35]
        /// </summary>
        private float _F0;

        /// <summary>
        /// // 垂直入射におけるフレネル反射係数の実部[フォトンマップ本 pp. 35]
        /// </summary>
        public float F0
        {
            get { return _F0; }
            set { _F0 = value; }
        }

        public TranslucentMaterial(float eta, float g, Vector3d sigma_s, Vector3d sigma_t, Vector3d rambda, float F0)
        {
            _eta = eta;
            _g = g;
            _sigma_s = sigma_s;
            _sigma_t = sigma_t;
            _rambda = rambda;
            _F0 = F0;
        }

        /// <summary>
        /// [Jensen2001]記載の数値から半透明物体の材質パラメータを設定する。
        /// </summary>
        /// <param name="sigma_s_dash">Reduced scattering coefficient. eg.(0.74, 0.88, 1.01) for Skin1</param>
        /// <param name="sigma_a">Absorption coefficient. eg.(0.032, 0.17, 0.48) for Skin1</param>
        /// <param name="eta">Relative index of refraction eg. 1.3 for Skin1</param>
        /// <param name="g">Phase function. A constant phase function results in isotropic scattering (g = 0).</param>
        static public TranslucentMaterial FromJensen2001(Vector3d sigma_s_dash, Vector3d sigma_a, float eta, float g)
        {

            Vector3d sigma_t_dash = sigma_s_dash + sigma_a;   // [Jensen 2001;pp2の左下]
            Vector3d sigma_s = 1.0f / (1.0f - g) * sigma_s_dash;   // [Jensen 2001;pp2の左下]     Scattering coefficient
            Vector3d sigma_t = sigma_a + sigma_s; // [Jensen2001; pp2の左上]              Extinction coefficient

            // scattering albedo
            // フォトンマップ本 pp.146
            Vector3d Rambda = sigma_s / sigma_t;

            // 垂直入射におけるフレネル反射係数の実部[フォトンマップ本 pp. 35]
            float F0 = ((1.0f - eta) / (1.0f + eta)) * ((1.0f - eta) / (1.0f + eta));

            return new TranslucentMaterial(eta, g, sigma_s, sigma_t, Rambda, F0);
        }

        /// <summary>
        /// フレネルリフレクタンス：どんくらい反射するか。（シュリックの近似）[フォトンマップ本 pp.35 式(2.30)]
        /// </summary>
        /// <param name="cos_theta">入射角θ1のコサイン</param>
        /// <returns></returns>
        public float Fr(float cos_theta1)
        {
            // フレネルリフレクタンス：どんくらい反射するか。（シュリックの近似）[フォトンマップ本 pp.35 式(2.30)]
            float c = 1.0f - cos_theta1;
            float Fr = F0 + (1.0f - F0) * c * c * c * c * c;  //シュリックの近似式

            return Fr;
        }

        /// <summary>
        /// フレネルトランスミッタンス（どんくらい中に入ってくか。）[フォトンマップ本 pp.35 式(2.32)のすぐ下に書いていてある。]
        /// </summary>
        /// <param name="cos_theta1">入射角θ1のコサイン</param>
        /// <returns></returns>
        public float Ft(float cos_theta1)
        {
            //フレネルトランスミッタンス（どんくらい中に入ってくか。）[フォトンマップ本 pp.35 式(2.32)のすぐ下に書いていてある。]
            float Ft = 1.0f- Fr(cos_theta1);
            return Ft;
        }


        public static TranslucentMaterial Apple { get { return FromJensen2001(new Vector3d(2.29, 2.39, 1.97), new Vector3d(0.003, 0.0034, 0.046), 1.3, 0.0); } }
        public static TranslucentMaterial Chicken1 { get { return FromJensen2001(new Vector3d(0.15, 0.21, 0.38), new Vector3d(0.015, 0.077, 0.19), 1.3, 0.0); } }
        public static TranslucentMaterial Chicken2 { get { return FromJensen2001(new Vector3d(0.19, 0.25, 0.32), new Vector3d(0.018, 0.088, 0.2), 1.3, 0.0); } }
        public static TranslucentMaterial Cream { get { return FromJensen2001(new Vector3d(7.38, 5.47, 3.15), new Vector3d(0.0002, 0.0028, 0.0163), 1.3, 0.0); } }
        public static TranslucentMaterial Ketchup { get { return FromJensen2001(new Vector3d(0.18, 0.07, 0.03), new Vector3d(0.061, 0.97, 1.45), 1.3, 0.0); } }
        public static TranslucentMaterial Marble { get { return FromJensen2001(new Vector3d(2.19, 2.62, 3), new Vector3d(0.0021, 0.0041, 0.0071), 1.5, 0.0); } }
        public static TranslucentMaterial Potato { get { return FromJensen2001(new Vector3d(0.68, 0.7, 0.55), new Vector3d(0.0024, 0.009, 0.12), 1.3, 0.0); } }
        public static TranslucentMaterial Skimmilk { get { return FromJensen2001(new Vector3d(0.7, 1.22, 1.9), new Vector3d(0.0014, 0.0025, 0.0142), 1.3, 0.0); } }
        public static TranslucentMaterial Skin1 { get { return FromJensen2001(new Vector3d(0.74, 0.88, 1.01), new Vector3d(0.032, 0.17, 0.48), 1.3, 0.0); } }
        public static TranslucentMaterial Skin2 { get { return FromJensen2001(new Vector3d(1.09, 1.59, 1.79), new Vector3d(0.013, 0.07, 0.145), 1.3, 0.0); } }
        public static TranslucentMaterial Spectralon { get { return FromJensen2001(new Vector3d(11.6, 20.4, 14.9), new Vector3d(0, 0, 0), 1.3, 0.0); } }
        public static TranslucentMaterial Wholemilk { get { return FromJensen2001(new Vector3d(2.55, 3.21, 3.77), new Vector3d(0.0011, 0.0024, 0.014), 1.3, 0.0); } }

    }
}
