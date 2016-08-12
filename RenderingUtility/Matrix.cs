//
//  Matrix.cs
//
//  Created by Hiroyuki Kubo on 2012/03/24.
//

using System;


namespace MathUtil
{
    public struct Matrix
    {
        private int _row;                   // 行の数
        private int _column;                // 列の数
        private float[,] _Element;                  //行列本体

        /// <summary>
        /// インスタンスコンストラクタ、行列のサイズが引数
        /// </summary>
        /// <param name="Row"></param>
        /// <param name="Column"></param>
        public Matrix(int Row, int Column)        
        {
            _row = Row;
            _column = Column;
            _Element = new float[_row, _column];

            for (int j = 0; j < _column; j++)       //全要素を0で初期化
            {
                for (int i = 0; i < _row; i++)
                {
                    _Element[i, j] = 0;
                }
            }
        }
        /// <summary>
        /// コピーコンストラクタ
        /// </summary>
        /// <param name="m"></param>
        public Matrix(Matrix m)
        {
            _row = m.Row;
            _column = m.Column;

            _Element = new float[_row, _column];

            for (int i = 0; i < _row; i++)
            {
                for (int j = 0; j < _column; j++)
                {
                    _Element[i, j] = m._Element[i, j];
                }
            }
        }

        /// <summary>
        /// 配列風に値を取り出す。（インデクサ）
        /// </summary>
        /// <param name="r"></param>
        /// <param name="c"></param>
        /// <returns></returns> 
        public float this[int r, int c]
        {
            get
            {
                return _Element[r, c];
            }
            set
            {
                _Element[r, c] = value;
            }
        }

        /// <summary>
        /// 行数プロパティ
        /// </summary>
        public int Row
        {
            set
            {
                _row = value;
            }
            get
            {
                return _row;
            }
        }

        /// <summary>
        /// 列数プロパティ
        /// </summary>
        public int Column
        {
            set
            {
                _column = value;
            }
            get
            {
                return _column;
            }
        }

        /// <summary>
        /// 行列同士の積
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        public static Matrix operator *(Matrix m1, Matrix m2)   //行列クラスの掛け算を＊で定義
        {
            Matrix result = new Matrix(m1._row, m2._column);

            if (m1._column == m2._row)
            {
                for (int i = 0; i < m1._row; i++)
                {
                    for (int j = 0; j < m2._column; j++)
                    {
                        for (int k = 0; k < m1._column; k++)
                        {
                            result._Element[i, j] += m1._Element[i, k] * m2._Element[k, j];
                        }
                    }
                }
            }
            else
            {
                throw new Exception("行列のサイズが合わず乗算できません。");
            }

            return result;
        }

        // 行列とスカラーの積
        public static Matrix operator *(Matrix m, float k)   //行列クラスの掛け算を＊で定義
        {
            int i, j;
            Matrix result = new Matrix(m._row, m._column);

            for (j = 0; j < m._column; j++) //列
            {
                for (i = 0; i < m._row; i++)  //行
                {
                    result._Element[i, j] = m._Element[i, j] * k;
                }
            }
            return result;
        }
        // 行列とスカラーの積
        public static Matrix operator *(float k, Matrix m)   //行列クラスの掛け算を＊で定義
        {
            int i, j;

            Matrix result = new Matrix(m._row, m._column);

            for (j = 0; j < m._column; j++) //列
            {
                for (i = 0; i < m._row; i++)  //行
                {
                    result._Element[i, j] = m._Element[i, j] * k;
                }
            }
            return result;
        }

        public static Matrix operator +(Matrix m1, Matrix m2)    //行列クラスの足し算を＋で定義
        {
            int i, j;
            Matrix result = new Matrix(m1._row, m1._column);

            if ((m1._row == m2._row) & (m1._column == m2._column))
            {
                for (j = 0; j < m1._column; j++)
                {
                    for (i = 0; i < m1._row; i++)
                    {
                        result._Element[i, j] = m1._Element[i, j] + m2._Element[i, j];
                    }
                }
            }
            else
            {
                throw new Exception("行列のサイズが合わず加算できません");
            }
            return result;
        }

        public static Matrix operator -(Matrix m1, Matrix m2)  //行列の引き算を-で定義
        {
            int i, j;
            Matrix result = new Matrix(m1._row, m1._column);

            if ((m1._row == m2._row) & (m1._column == m2._column))
            {
                for (j = 0; j < m1._column; j++)
                {
                    for (i = 0; i < m1._row; i++)
                    {
                        result._Element[i, j] = m1._Element[i, j] - m2._Element[i, j];
                    }
                }
            }
            else
            {
                throw new Exception("行列のサイズが合わず減算できません");
            }
            return result;
        }

        public static Matrix CopyFrom(Matrix m)
        {
            return new Matrix(m);
        }

        /// <summary>
        /// 転置行列計算メソッド
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static Matrix Transpose(Matrix m)
        {
            int i, j;
            Matrix result = new Matrix(m._column, m._row);
            for (j = 0; j < m._column; j++)
            {
                for (i = 0; i < m._row; i++)
                {
                    result._Element[j, i] = m._Element[i, j];
                }
            }
            return result;
        }

        public override string ToString()
        {
            string res = "[" + Row + "," + _column + "]";
            return res;
        }

        ///http://lauelab.unh.edu/projects/AOS/browser/trunk/Numerical%20recipes/ansi_c/ludcmp.c?rev=2438
        private static void ludcmp(float[,] a, int n, int[] indx, ref float d)
        {
            int i, imax, j, k;
            float big, dum, sum, temp;
            float[] vv = new float[n];						// 各行の暗黙のスケーリングを記録する．
            const float TINY = 1.0e-20f;
            imax = -10;

            d = 1.0f;							// まだ行交換していない．
            for (i = 0; i < n; i++)
            {				// 行についてループし，暗黙のスケーリングの情報を得る．
                big = 0.0f;
                for (j = 0; j < n; j++)
                    if ((temp = Math.Abs(a[i, j])) > big) big = temp;
                if (big == 0.0) throw new Exception("Singular matrix in routine ludcmp\n");	// 最大要素が０なら特異行列である．
                vv[i] = 1.0f / big;				// スケーリングを記録する．
            }
            for (j = 0; j < n; j++)
            {
                // Crout法，列についてのループ
                for (i = 0; i < j; i++)
                {
                    // 方程式(2.3.12)のi=j以外
                    sum = a[i, j];
                    for (k = 0; k < i; k++) sum -= a[i, k] * a[k, j];
                    a[i, j] = sum;
                }
                big = 0.0f;
                for (i = j; i < n; i++)
                {
                    sum = a[i, j];
                    for (k = 0; k < j; k++)
                        sum -= a[i, k] * a[k, j];
                    a[i, j] = sum;
                    if ((dum = vv[i] * Math.Abs(sum)) >= big)
                    {
                        big = dum;
                        imax = i;
                    }
                }
                if (j != imax)
                {
                    for (k = 0; k < n; k++)
                    {
                        dum = a[imax, k];
                        a[imax, k] = a[j, k];
                        a[j, k] = dum;
                    }
                    d = -d;
                    vv[imax] = vv[j];
                }
                indx[j] = imax;
                if (a[j, j] == 0.0) a[j, j] = TINY;
                if (j != n)
                {
                    dum = 1.0f / (a[j, j]);
                    for (i = j + 1; i < n; i++) a[i, j] *= dum;
                }
            }
            //free_vector(vv);
        }
        private static void lubksb(float[,] a, int n, int[] indx, float[] b)
        {
            int i, ii = 0, ip, j;
            float sum;

            for (i = 0; i < n; i++)
            {
                ip = indx[i];
                sum = b[ip];
                b[ip] = b[i];
                if (ii != 0)
                    for (j = ii - 1; j <= i - 1; j++) sum -= a[i, j] * b[j];
                else if (sum != 0.0) ii = i + 1;
                b[i] = sum;
            }
            for (i = n - 1; i >= 0; i--)
            {
                sum = b[i];
                for (j = i + 1; j < n; j++) sum -= a[i, j] * b[j];
                b[i] = sum / a[i, i];
            }
        }
        /// <summary>
        /// 逆行列を求める。（引数は破壊されない。）
        /// </summary>
        /// <param name="m1"></param>
        /// <returns></returns>
        public static Matrix Inverse(Matrix m1)
        {
            Matrix result = new Matrix(m1._row, m1._column);

            int N = m1._row;
            int[] indx = new int[N];
            float d = 0;
            float[] col = new float[N];

            float[,] data = new float[N, N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    data[i, j] = m1._Element[i, j];
                }
            }

            ludcmp(data, N, indx, ref d);

            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++) col[i] = 0.0f;
                col[j] = 1.0f;
                lubksb(data, N, indx, col);

                for (int i = 0; i < N; i++) result._Element[i, j] = col[i];

            }
            return result;
        }

        /// <summary>
        /// 疑似逆行列を求める（引数は破壊されない。）
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static Matrix PesudoInverse(Matrix m)
        {
            Matrix U, V, S;
            SVD(m, out U, out V, out S);

            return V * Matrix.Inverse(S) * Matrix.Transpose(U);
        }

        /// <summary>
        /// M = U * S * V^t に分解するルーチン。
        /// S:実数対角行列
        /// U,V：直交またはユニタリ正方行列
        /// </summary>
        /// <param name="M"></param>
        /// <param name="U"></param>
        /// <param name="V"></param>
        /// <param name="S"></param>
        public static void SVD(Matrix M, out Matrix U, out Matrix V, out Matrix S)
        {
            int m = M.Row;
            int n = M.Column;

            float[,] a = new float[m, n];

            float[] s = new float[n]; // 対角成分
            float[,] v = new float[n, n];

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i, j] = M[i, j];
                }
            }


            svdcmp(a, m, n, s, v);

            U = new Matrix(m, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    U[i, j] = a[i, j];
                }
            }

            V = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    V[i, j] = v[i, j];
                }
            }

            S = new Matrix(n, n);
            for (int i = 0; i < n; i++)
            {
                S[i, i] = s[i];

            }
        }

        /// <summary>
        /// http://www-bs.ss.oka-pu.ac.jp/inoue/numerical/svdcmp_poly_fit2.txt
        /// A = U * W * V^t に分解するルーチン。Uはaに上書きされる。
        /// </summary>
        /// <param name="a"></param>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="w"></param>
        /// <param name="v"></param>
        static void svdcmp(float[,] a, int m, int n, float[] w, float[,] v)
        {

            int flag, i, its, j, jj, k, l, nm;
            float anorm, c, f, g, h, s, scale, x, y, z;
            float[] rv1 = new float[n];
            l = nm = 0;


            g = scale = anorm = 0.0f;
            for (i = 0; i < n; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = s = scale = 0.0f;
                if (i < m)
                {
                    for (k = i; k < m; k++) scale += Math.Abs(a[k, i]);
                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            a[k, i] /= scale;
                            s += a[k, i] * a[k, i];
                        }
                        f = a[i, i];
                        g = -SIGN((float)Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, i] = f - g;
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0f, k = i; k < m; k++) s += a[k, i] * a[k, j];
                            f = s / h;
                            for (k = i; k < m; k++) a[k, j] += f * a[k, i];
                        }
                        for (k = i; k < m; k++) a[k, i] *= scale;
                    }
                }
                w[i] = scale * g;
                g = s = scale = 0.0f;
                if (i < m && i != n)
                {
                    for (k = l; k < n; k++) scale += Math.Abs(a[i, k]);
                    if (scale != 0.0)
                    {
                        for (k = l; k < n; k++)
                        {
                            a[i, k] /= scale;
                            s += a[i, k] * a[i, k];
                        }
                        f = a[i, l];
                        g = -SIGN((float)Math.Sqrt(s), f);
                        h = f * g - s;
                        a[i, l] = f - g;
                        for (k = l; k < n; k++) rv1[k] = a[i, k] / h;
                        for (j = l; j < m; j++)
                        {
                            for (s = 0.0f, k = l; k < n; k++) s += a[j, k] * a[i, k];
                            for (k = l; k < n; k++) a[j, k] += s * rv1[k];
                        }
                        for (k = l; k < n; k++) a[i, k] *= scale;
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(w[i]) + Math.Abs(rv1[i])));
            }
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                            v[j, i] = (a[i, j] / a[i, l]) / g;
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0f, k = l; k < n; k++) s += a[i, k] * v[k, j];
                            for (k = l; k < n; k++) v[k, j] += s * v[k, i];
                        }
                    }
                    for (j = l; j < n; j++) v[i, j] = v[j, i] = 0.0f;
                }
                v[i, i] = 1.0f;
                g = rv1[i];
                l = i;
            }
            for (i = (int)Math.Min(m, n) - 1; i >= 0; i--)
            {
                l = i + 1;
                g = w[i];
                for (j = l; j < n; j++) a[i, j] = 0.0f;
                if (g != 0.0)
                {
                    g = 1.0f / g;
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0f, k = l; k < m; k++) s += a[k, i] * a[k, j];
                        f = (s / a[i, i]) * g;
                        for (k = i; k < m; k++) a[k, j] += f * a[k, i];
                    }
                    for (j = i; j < m; j++) a[j, i] *= g;
                }
                else for (j = i; j < m; j++) a[j, i] = 0.0f;
                ++a[i, i];
            }
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 1; its <= 30; its++)
                {
                    flag = 1;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l;
                        if ((float)(Math.Abs(rv1[l]) + anorm) == anorm)
                        {
                            flag = 0;
                            break;
                        }
                        if ((float)(Math.Abs(w[nm]) + anorm) == anorm) break;
                    }
                    if (flag != 0)
                    {
                        c = 0.0f;
                        s = 1.0f;
                        for (i = l; i <= k; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if ((float)(Math.Abs(f) + anorm) == anorm) break;
                            g = w[i];
                            h = pythag(f, g);
                            w[i] = h;
                            h = 1.0f / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < m; j++)
                            {
                                y = a[j, nm];
                                z = a[j, i];
                                a[j, nm] = y * c + z * s;
                                a[j, i] = z * c - y * s;
                            }
                        }
                    }
                    z = w[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k] = -z;
                            for (j = 0; j < n; j++) v[j, k] = -v[j, k];
                        }
                        break;
                    }
                    if (its == 30) throw new Exception("no convergence in 30 svdcmp iterations");
                    x = w[l];
                    nm = k;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
                    g = pythag(f, 1.0f);
                    f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                    c = s = 1.0f;
                    for (j = l; j < nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = pythag(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 0; jj < n; jj++)
                        {
                            x = v[jj, j];
                            z = v[jj, i];
                            v[jj, j] = x * c + z * s;
                            v[jj, i] = z * c - x * s;
                        }
                        z = pythag(f, h);
                        w[j] = z;
                        if (z != 0.0)
                        {
                            z = 1.0f / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 0; jj < m; jj++)
                        {
                            y = a[jj, j];
                            z = a[jj, i];
                            a[jj, j] = y * c + z * s;
                            a[jj, i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0f;
                    rv1[k] = f;
                    w[k] = x;
                }
            }
        }

        public static float Determinant(Matrix m)
        {
            if (m.Row != m.Column) throw new Exception("行列式は正方行列でないと求まらない。");

            Matrix result = new Matrix(m._row, m._column);

            int N = m._row;
            float d = 0;
            int[] indx = new int[N];

            float[,] data = new float[N, N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    data[i, j] = m._Element[i, j];
                }
            }

            ludcmp(data, N, indx, ref d);

            for (int i = 0; i < N; i++) d *= data[i, i];

            return d;
        }

        private static float pythag(float a, float b)
        {
            float absa, absb;
            absa = Math.Abs(a);
            absb = Math.Abs(b);
            if (absa > absb) return absa * (float)Math.Sqrt(1.0 + SQR(absb / absa));
            else return (absb == 0.0f ? 0.0f : absb * (float)Math.Sqrt(1.0 + SQR(absa / absb)));
        }
        private static float SQR(float val) { return val * val; }
        private static float SIGN(float a, float b)
        {
            return b >= 0.0 ? Math.Abs(a) : -Math.Abs(a);
        }


        /// <summary>
        /// 実数対象行列の固有値と固有ベクトルを求める関数
        /// </summary>
        /// <param name="mat">入力となる実対象行列（破壊されない）</param>
        /// <param name="eval">出力となる固有値</param>
        /// <param name="evec">対応する固有ベクトル</param>
        public static void EigenRS(Matrix mat, out float[] eval, out float[][] evec)
        {
            if(mat.Row != mat.Column) throw new Exception("正方行列を入力して下さい。");

            // 行列のサイズ
            int n = mat.Row;

            // まずは値をコピー
            float[,] val = new float[n + 1, n + 1];   // Numerical Recipesの実装なので[0]要素は使わない。

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    val[i + 1, j + 1] = mat[i, j];
                }
            }

            float[] e = new float[n + 1];
            float[] d = new float[n + 1];

            // 三重対角化
            tred2(val, n, d, e);

            // QL法
            tqli(d, e, n, val);

            // 固有値・固有ベクトルのコピー
            eval = new float[n];
            for (int i = 0; i < n; i++) eval[i] = d[i + 1];

            evec = new float[n][];
            for (int i = 0; i < n; i++)
            {
                evec[i] = new float[n];
                for (int j = 0; j < n; j++) evec[i][j] = val[j + 1, i + 1];
            }

            // 固有値が降順に並び替え
            Array.Sort(eval, evec); // まず昇順
            Array.Reverse(eval);    // 反転させて降順にする
            Array.Reverse(evec);
        }

        /// <summary>
        /// 単位行列を生成する
        /// </summary>
        /// <param name="d">行列のサイズ（d次の正方行列の単位行列が生成される。）</param>
        /// <returns></returns>
        public static Matrix Identity(int d)
        {
            Matrix m = new Matrix(d, d);

            for (int i = 0; i < d; i++) m[d, d] = 1.0f;

            return m;
        }

        private static void tqli(float[] d, float[] e, int n, float[,] z)
        {
            //float pythag(float a, float b);
            int m, l, iter, i, k;
            float s, r, p, g, f, dd, c, b;
            for (i = 2; i <= n; i++) e[i - 1] = e[i]; //Convenient to renumber the elements of e.
            e[n] = 0.0f;
            for (l = 1; l <= n; l++)
            {
                iter = 0;
                do
                {
                    for (m = l; m <= n - 1; m++)
                    { //Look for a single small subdiagonal element to split the matrix.
                        dd = Math.Abs(d[m]) + Math.Abs(d[m + 1]);
                        if ((float)(Math.Abs(e[m]) + dd) == dd) break;
                    }
                    if (m != l)
                    {
                        if (iter++ == 30) throw new Exception("Too many iterations in tqli");
                        g = (d[l + 1] - d[l]) / (2.0f * e[l]); //Form shift.
                        r = pythag(g, 1.0f);
                        g = d[m] - d[l] + e[l] / (g + SIGN(r, g)); //This is dm − ks.
                        s = c = 1.0f;
                        p = 0.0f;
                        for (i = m - 1; i >= l; i--)
                        { //A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
                            f = s * e[i];
                            b = c * e[i];
                            e[i + 1] = (r = pythag(f, g));
                            if (r == 0.0)
                            { //Recover from underflow.
                                d[i + 1] -= p;
                                e[m] = 0.0f;
                                break;
                            }
                            s = f / r;
                            c = g / r;
                            g = d[i + 1] - p;
                            r = (d[i] - g) * s + 2.0f * c * b;
                            d[i + 1] = g + (p = s * r);
                            g = c * r - b;
                            /* Next loop can be omitted if eigenvectors not wanted*/
                            for (k = 1; k <= n; k++)
                            { //Form eigenvectors.
                                f = z[k, i + 1];
                                z[k, i + 1] = s * z[k, i] + c * f;
                                z[k, i] = c * z[k, i] - s * f;
                            }
                        }
                        if (r == 0.0 && i >= l) continue;
                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0f;
                    }
                } while (m != l);
            }
        }
        private static void tred2(float[,] a, int n, float[] d, float[] e)
        {
            int l, k, j, i;
            float scale, hh, h, g, f;
            for (i = n; i >= 2; i--)
            {
                l = i - 1;
                h = scale = 0.0f;
                if (l > 1)
                {
                    for (k = 1; k <= l; k++)
                        scale += Math.Abs(a[i,k]);
                    if (scale == 0.0) //Skip transformation.
                        e[i] = a[i,l];
                    else
                    {
                        for (k = 1; k <= l; k++)
                        {
                            a[i,k] /= scale; //Use scaled a’s for transformation.
                            h += a[i,k] * a[i,k];// Form σ in h.
                        }
                        f = a[i,l];
                        g = (f >= 0.0f ? -(float)Math.Sqrt(h) : (float)Math.Sqrt(h));
                        e[i] = scale * g;
                        h -= f * g; //Now h is equation (11.2.4).
                        a[i,l] = f - g; //Store u in the ith row of a.
                        f = 0.0f;
                        for (j = 1; j <= l; j++)
                        {
                            /* Next statement can be omitted if eigenvectors not wanted */
                            a[j,i] = a[i,j] / h; //Store u/H in ith column of a.
                            g = 0.0f; //Form an element of A · u in g.
                            for (k = 1; k <= j; k++)
                                g += a[j,k] * a[i,k];
                            for (k = j + 1; k <= l; k++)
                                g += a[k,j] * a[i,k];
                            e[j] = g / h; //Form element of p in temporarily unused element of e.
                            f += e[j] * a[i,j];
                        }
                        hh = f / (h + h); //Form K, equation (11.2.11).
                        for (j = 1; j <= l; j++)
                        { //Form q and store in e overwriting p.
                            f = a[i,j];
                            e[j] = g = e[j] - hh * f;
                            for (k = 1; k <= j; k++) //Reduce a, equation (11.2.13).
                                a[j,k] -= (f * e[k] + g * a[i,k]);
                        }
                    }
                }
                else
                    e[i] = a[i,l];
                d[i] = h;
            }
            /* Next statement can be omitted if eigenvectors not wanted */
            d[1] = 0.0f;
            e[1] = 0.0f;
            /* Contents of this loop can be omitted if eigenvectors not
            wanted except for statement d[i]=a[i,i]; */
            for (i = 1; i <= n; i++)
            { //Begin accumulation of transformation matrices.
                l = i - 1;
                if (d[i] != 0.0)
                { //This block skipped when i=1.
                    for (j = 1; j <= l; j++)
                    {
                        g = 0.0f;
                        for (k = 1; k <= l; k++) //Use u and u/H stored in a to form P·Q.
                            g += a[i,k] * a[k,j];
                        for (k = 1; k <= l; k++)
                            a[k,j] -= g * a[k,i];
                    }
                }
                d[i] = a[i,i]; //This statement remains.
                a[i,i] = 1.0f; //Reset row and column of a to identity
                for (j = 1; j <= l; j++) a[j,i] = a[i,j] = 0.0f; // for next iteration.
            }
        }
    }
}
    
