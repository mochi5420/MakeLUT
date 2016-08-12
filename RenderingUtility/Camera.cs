using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class Camera
    {

        private Vector3f _xaxis, _yaxis, _zaxis;
        private float _aspect, _vfov;
        private float _near, _far;
        private Vector3f _pos;


        public Camera(Vector3f pos, Vector3f target, Vector3f up, float aspect, float vfov, float nearclip, float farclip)
        {
            // カメラから見た奥行き方向
            _zaxis = Vector3f.Normalize(target - pos);

            // カメラから見た上向き方向
            _yaxis = Vector3f.Normalize(up - (Vector3f.Dot(up, _zaxis) * _zaxis));    // upベクトルが奥行き方向と直交しないときの対処

            // 上記１つの方向に直交する方向
            _xaxis = Vector3f.Normalize(Vector3f.Cross(_zaxis, _yaxis));


            _pos = pos;
            _aspect = aspect;
            _vfov = vfov;
            _near = nearclip;
            _far = farclip;

        }

        /// <summary>
        /// 注目するピクセルをレンダリングするためのレイを計算する
        /// </summary>
        /// <param name="x">画面座標系の横軸成分[-1, +1]</param>
        /// <param name="y">画面座標系の縦軸成分[-1, +1]</param>
        /// <returns></returns>
        public Ray CameraRay(float x, float y)
        {
            // 注目するピクセルの視線の方向
            Vector3f r = _zaxis + (float)Math.Tan(0.5 * _vfov) * x * _xaxis + 1.0f / _aspect * (float)Math.Tan(0.5 * _vfov) * y * _yaxis;

            // 
            Ray ray = new Ray(_pos + _near * r, (_far - _near) * r);

            return ray;
        }
    }
}
