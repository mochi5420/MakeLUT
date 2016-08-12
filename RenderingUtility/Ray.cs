using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class Ray
    {
        private Vector3f _center;
        private Vector3f _half;
        private Vector3f _dir;  // 正規化されたDirection

        public Ray(Vector3f Origin, Vector3f Direction)
        {
            _center = Origin + 0.5f * Direction;
            _half = 0.5f * Direction;
            _dir = Vector3f.Normalize(Direction);
        }
        public void FromOriginDirection(Vector3f Origin, Vector3f Direction)
        {
            _center = Origin + 0.5f * Direction;
            _half = 0.5f * Direction;
            _dir = Vector3f.Normalize(Direction);
        }
        public static Ray FromPoints(Vector3f A, Vector3f B)
        {
            Vector3f origin = A;
            Vector3f dir = B - A;
            return new Ray(origin, dir);
        }
        public Vector3f Center
        {
            get { return _center; }
            set { _center = value; }
        }
        public Vector3f HalfVector
        {
            get { return _half; }
            set { _half = value; }
        }
        // 正規化済みの方向ベクトル
        public Vector3f Direction
        {
            get { return _dir; }
        }
        public Vector3f Origin
        {
            get { return _center - _half; }
        }
        /// <summary>
        /// 始点
        /// </summary>
        public Vector3f PosA
        {
            get { return _center - _half; }
        }
        /// <summary>
        /// 終点
        /// </summary>
        public Vector3f PosB
        {
            get { return _center + _half; }
        }

    }
}
