using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{
    public class Sphere
    {
        private float _radius;

        public float Radius
        {
            get { return _radius; }
            set { _radius = value; }
        }
        private Vector3f _center;

        public Vector3f Center
        {
            get { return _center; }
            set { _center = value; }
        }


        public Sphere(float r, Vector3f center)
        {
            _radius = r;
            _center = center;
        }
    }
}