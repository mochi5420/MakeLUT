using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using Utility;

namespace RenderingUtility
{
    public class Texture2D
    {

        private float[, ,] _imagedata;
        private int _width;
        private int _height;

        public static Texture2D FromBitmap(Bitmap bmp)
        {
            int height = bmp.Height;
            int width = bmp.Width;

            float[, ,] data = new float[height, width, 3];

            for (int h = 0; h < height; h++)
            {
                for (int w = 0; w < width; w++)
                {
                    Color c = bmp.GetPixel(w, h);
                    data[h, w, 0] = (float)c.R / 255f;
                    data[h, w, 1] = (float)c.G / 255f;
                    data[h, w, 2] = (float)c.B / 255f;
                }
            }

            return new Texture2D(data);
        }

        public Texture2D(float[, ,] image)
        {
            _imagedata = (float[,,])image.Clone();
            _height = _imagedata.GetLength(0);
            _width = _imagedata.GetLength(1);

        }

        private float SampleOne(int w, int h, int ch)
        {
            int w2 = w % _width;
            int h2 = h % _height;

            if (w2 < 0) w2 += _width;
            if (h2 < 0) h2 += _width;

            return _imagedata[h2, w2, ch];
        }

        /// <summary>
        /// http://imagingsolution.blog107.fc2.com/blog-entry-142.html
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <param name="ch"></param>
        /// <returns></returns>
        public float SampleOne(float u, float v, int ch)
        {
            float x = u * (float)_width;
            float y = v * (float)_height;

            // fx ... [x]
            float fx = (float)Math.Floor(x);
            float fy = (float)Math.Floor(y);

            return (fx + 1f - x) * (fy + 1f - y) * SampleOne((int)fx, (int)fy, ch)
                + (fx + 1f - x) * (y - fy) * SampleOne((int)(fx + 1f), (int)fy, ch)
                + (x - fx) * (fy + 1f - y) * SampleOne((int)fx, (int)(fy + 1f), ch)
                + (x - fx) * (y - fy) * SampleOne((int)(fx + 1f), (int)(fy + 1f), ch);
        }
        public Vector3f Sample(Vector2f tc)
        {
            return new Vector3f(this.SampleOne(tc.X, tc.Y, 0), this.SampleOne(tc.X, tc.Y, 1), this.SampleOne(tc.X, tc.Y, 2));
        }

    }
}
