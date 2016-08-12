using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Utility;
using RenderingUtility;
using MeshUtil;
using MathUtil;



namespace Raytrace_practice
{
    public partial class Form1 : Form
    {
        int sumpleNum = 500;
        WavefrontObj obj;
        Scene scene;

        public Form1()
        {
            InitializeComponent();
        }

        private static void WriteCsv(Vector3f[] Position)
        {
            try
            {
                using (var sw = new System.IO.StreamWriter(@"test.csv"))
                {
                    for (int i = 0; i < Position.Length; i++)
                    {
                        sw.WriteLine("{0}, {1}, {2}", Position[i].X, Position[i].Y, Position[i].Z);

                    }

                    //閉じる
                    sw.Close();
                    System.Console.WriteLine("Finish!");
                }


            }
            catch (System.Exception e)
            {
                // ファイルを開くのに失敗したときエラーメッセージを表示
                System.Console.WriteLine(e.Message);

            }

        }


        private float get_thickness2(Ray ray, Vector3f p0, Vector3f light_dir, float thickness)
        {

            float u, v, t;

            Triangle tri;

            if (scene.RayTrace(ray, out tri, out u, out v, out t) == IntersectionStatus.Intersect)
            {

                //衝突した点の法線
                Vector3f n = Vector3f.Normalize(tri.FaceNormal);

                //衝突した点の座標
                Vector3f p = ray.Origin + t * ray.Direction;

                //法線と進行方向の内積
                float Dot = Vector3f.Dot(n, light_dir);

                //厚み
                if (Dot < 0.0f)
                {
                    thickness = Vector3f.Length(p - p0);
                }
                else
                {
                    thickness = 0.0f;
                }

                //HITした点から再びRayを飛ばそう
                //面法線とRay方向が同じ方向なら＋、逆向きならー側に少しずらす
                float eps = 0.00001f;
                Vector3f p_dash;

                if (Dot < 0.0f)
                {
                    p_dash.X = p.X - n.X * eps;
                    p_dash.Y = p.Y - n.Y * eps;
                    p_dash.Z = p.Z - n.Z * eps;
                }
                else
                {
                    p_dash.X = p.X + n.X * eps;
                    p_dash.Y = p.Y + n.Y * eps;
                    p_dash.Z = p.Z + n.Z * eps;
                }

                //再びRay
                Ray ray2 = new Ray(p_dash, light_dir);

                thickness = thickness + get_thickness2(ray2, p, light_dir, thickness);
            }
            else thickness = 0.0f;

            return thickness;
        }

        private void button1_Click(object sender, EventArgs e)
        {

            //objファイルの読込（三角形メッシュのみ可）
            obj = new WavefrontObj("../../winebottle0.obj");

            //高速化のためにシーン作成
            scene = new Scene(new AABB(new Vector3f(0.009f, 0.009f, 0.009f), new Vector3f(0.5f, 0.5f, 0.5f)),
                                40, 40, 40);

            //objのための三角形メッシュを用意
            Triangle[] triangles = new Triangle[obj.NumFacet];

            //sceneに三角形メッシュ登録
            for (int i = 0; i < obj.NumFacet; i++)
            {
                Vector3d p0 = obj.Position(i, 0);
                Vector3d p1 = obj.Position(i, 1);
                Vector3d p2 = obj.Position(i, 2);

                triangles[i] = new Triangle(new Vector3f((float)p0.X, (float)p0.Y, (float)p0.Z),
                                                new Vector3f((float)p1.X, (float)p1.Y, (float)p1.Z),
                                                new Vector3f((float)p2.X, (float)p2.Y, (float)p2.Z));
    
                scene.RegistTriangle(triangles[i]);
            }

            //とりあえず適当なサンプリング。どうせあとでポアソンディスクサンプリングとかするから許して
            int interval = obj.NumFacet / sumpleNum;
            Vector3f[] samp_pos = new Vector3f[sumpleNum];

            for (int i = 0; i<sumpleNum; i++)
            {
                //i番目の三角形メッシュの真ん中あたり（クソ適当）
                samp_pos[i] = new Vector3f((float)((obj.Position(i * interval, 0).X + obj.Position(i * interval, 1).X + obj.Position(i * interval, 2).X) / 3.0),
                                                    (float)((obj.Position(i * interval, 0).Y + obj.Position(i * interval, 1).Y + obj.Position(i * interval, 2).Y) / 3.0),
                                                    (float)((obj.Position(i * interval, 0).Z + obj.Position(i * interval, 1).Z + obj.Position(i * interval, 2).Z) / 3.0));
                
            }

            
            //頂点の座標を書き出す
            WriteCsv(samp_pos);

            //縦横の分割数
            int width = pictureBox1.Width;
            int height = pictureBox1.Height;

            ////結果を格納するBitmap
            Bitmap _bitmap = null;
            try
            {
                _bitmap = new Bitmap(height, width);
            }
            catch (System.Exception ex)
            {
                Console.WriteLine(ex.Message);
            }

            //厚みを格納する配列の初期化
            float[,] thickness = new float[height, width];
            for (int h = 0; h < height; h++)
            {
                for (int w = 0; w < width; w++)
                {
                    thickness[h, w] = 0.0f;
                }
            }

            //厚みの計算。
            int count = 0;
            float progress = 0.0f;
            Parallel.For(0, height * width, i =>
            {
                int w = i % width;
                int h = ((i - i % width) / width) % height;

                if (h != w)
                {
                    if (thickness[h, w] == 0.0f)
                    {
                        Vector3f dir = samp_pos[w] - samp_pos[h];
                        dir = Vector3f.Normalize(dir);
                        float length = Vector3f.Length(samp_pos[w] - samp_pos[h]);

                        Ray ray = new Ray(samp_pos[h], dir);

                        //空洞(cavity)部の厚み取得
                        float cavity = 0.0f;
                        cavity = get_thickness2(ray, samp_pos[h], dir, cavity);

                        if (cavity == 0.0f)
                        {
                            thickness[h, w] = length;
                        }
                        else
                        {
                            thickness[h, w] = length - cavity;

                            if (thickness[h, w] < 0.0f)
                            {
                                thickness[h, w] = 0.0f;
                            }
                        }

                        thickness[w, h] = thickness[h, w];
                    }

                }

                else
                {
                    thickness[h, w] = 0.0f;
                }
                count += 1;
                progress = (float)count / (float)(sumpleNum * sumpleNum) * 100.0f;
                Console.WriteLine("{0}%" ,progress);

            });

            //厚みの最大値
            float t_max = float.MinValue;

            for (int h = 0; h < height; h++)
            {
                for (int w = 0; w < width; w++)
                {
                    if (t_max < thickness[h, w])
                    {
                        t_max = thickness[h, w];
                    }
                }
            }


            //輝度値を計算
            for (int h = 0; h < height; h++)
            {
                for (int w = 0; w < width; w++)
                {

                    //厚みを０～２５５のカラーに収める
                    float _color = 255.0f * thickness[h, w] / t_max;

                    // bmpに格納
                    _bitmap.SetPixel(w, h, Color.FromArgb((int)_color, (int)_color, (int)_color));
                    //_bitmap.SetPixel(w, height - h - 1, Color.FromArgb((int)_color, (int)_color, (int)_color));
                }
            }

            //pictureboxに渡す
            pictureBox1.Image = _bitmap;

            pictureBox1.Image.Save("winebottle0_test.png");
        }
    }
}


