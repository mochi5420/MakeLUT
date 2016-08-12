using System;
using System.Collections.Generic;
using System.Text;
using Utility;

namespace RenderingUtility
{

    public enum IntersectionStatus
    {
        Reject,
        Intersect,
    };
    public enum OverlapStatus
    {
        Disjoint,
        Overlap,
    };

    public static class Intersections
    {
        private static float epsilon = float.Epsilon;

        private static bool AXISTEST_X01(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {
            float min, max;
            float p0 = a * v0.Y - b * v0.Z;
            float p2 = a * v2.Y - b * v2.Z;

            if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }

            float rad = fa * boxhalfsize.Y + fb * boxhalfsize.Z;

            if (min > rad || max < -rad) return false;
            else return true;
        }
        private static bool AXISTEST_X2(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {

            float min, max;
            float p0 = a * v0.Y - b * v0.Z;
            float p1 = a * v1.Y - b * v1.Z;

            if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }

            float rad = fa * boxhalfsize.Y + fb * boxhalfsize.Z;

            if (min > rad || max < -rad) return false;
            else return true;
        }
        private static bool AXISTEST_Y02(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {
            float min, max;
            float p0 = -a * v0.X + b * v0.Z;
            float p2 = -a * v2.X + b * v2.Z;

            if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }

            float rad = fa * boxhalfsize.X + fb * boxhalfsize.Z;

            if (min > rad || max < -rad) return false;
            else return true;

        }
        private static bool AXISTEST_Y1(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {
            float min, max;
            float p0 = -a * v0.X + b * v0.Z;
            float p1 = -a * v1.X + b * v1.Z;

            if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }

            float rad = fa * boxhalfsize.X + fb * boxhalfsize.Z;

            if (min > rad || max < -rad) return false;
            else return true;

        }
        private static bool AXISTEST_Z12(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {
            float min, max;
            float p1 = a * v1.X - b * v1.Y;
            float p2 = a * v2.X - b * v2.Y;

            if (p2 < p1) { min = p2; max = p1; } else { min = p1; max = p2; }

            float rad = fa * boxhalfsize.X + fb * boxhalfsize.Y;

            if (min > rad || max < -rad) return false;
            else return true;

        }
        private static bool AXISTEST_Z0(float a, float b, float fa, float fb, Vector3f v0, Vector3f v1, Vector3f v2, Vector3f boxhalfsize)
        {
            float min, max;
            float p0 = a * v0.X - b * v0.Y;
            float p1 = a * v1.X - b * v1.Y;

            if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }

            float rad = fa * boxhalfsize.X + fb * boxhalfsize.Y;

            if (min > rad || max < -rad) return false;
            else return true;
        }
        private static void FINDMINMAX(float x0, float x1, float x2, out float min, out float max)
        {
            min = max = x0;
            if (x1 < min) min = x1;
            if (x1 > max) max = x1;
            if (x2 < min) min = x2;
            if (x2 > max) max = x2;
        }

        /// <summary>
        /// RayとShereの衝突判定（リアルタイムレンダリング第2版・日本語版・pp.487）
        /// </summary>
        /// <param name="ray">レイ</param>
        /// <param name="s">球</param>
        /// <param name="t">レイの原点から交点までの距離</param>
        /// <param name="p">交点（衝突する場合）</param>        /// <returns></returns>
        public static IntersectionStatus RaySphereIntersect(Ray ray, Sphere s, out float t, out Vector3f p)
        {
            return RaySphereIntersect(ray.Origin, ray.Direction, s.Center, s.Radius, out t, out p);
        }
        /// <summary>
        /// RayとShereの衝突判定（リアルタイムレンダリング第2版・日本語版・pp.487）
        /// </summary>
        /// <param name="ray">レイ</param>
        /// <param name="s">球</param>
        /// <param name="p">交点（衝突する場合）</param>        /// <returns></returns>
        public static IntersectionStatus RaySphereIntersect(Ray ray, Sphere s, out Vector3f p)
        {
            float t;
            return RaySphereIntersect(ray.Origin, ray.Direction, s.Center, s.Radius, out t, out p);
        }
        /// <summary>
        /// RayとShereの衝突判定（リアルタイムレンダリング第2版・日本語版・pp.487）
        /// </summary>
        /// <param name="ray">レイ</param>
        /// <param name="c">球の中心</param>
        /// <param name="r">球の半径</param>
        /// <param name="t">レイの原点から交点までの距離</param>
        /// <param name="p">交点（衝突する場合）</param>
        /// <returns></returns>
        public static IntersectionStatus RaySphereIntersect(Ray ray, Vector3f c, float r, out float t, out Vector3f p)
        {
            return RaySphereIntersect(ray.Origin, ray.Direction, c, r, out t, out p);
        }
        /// <summary>
        /// RayとShereの衝突判定（リアルタイムレンダリング第2版・日本語版・pp.487）
        /// </summary>
        /// <param name="o">レイの原点</param>
        /// <param name="d">レイの方向</param>
        /// <param name="c">球の中心</param>
        /// <param name="r">球の半径</param>
        /// <param name="t">レイの原点から交点までの距離</param>
        /// <param name="p">交点（衝突する場合）</param>
        /// <returns></returns>
        public static IntersectionStatus RaySphereIntersect(Vector3f o, Vector3f d, Vector3f c, float r,out float t, out Vector3f p)
        {
            Vector3f l = c - o;
            float s = Vector3f.Dot(l, d);
            float l_sq = Vector3f.Dot(l, l);
            if (s < 0.0 && l_sq > r * r)
            {
                t = 0;
                p = new Vector3f();
                return IntersectionStatus.Reject;
            }

            float m_sq = l_sq - s * s;
            if (m_sq > r * r)
            {
                t = 0;
                p = new Vector3f();
                return IntersectionStatus.Reject;
            }

            float q = (float)Math.Sqrt(r * r - m_sq);
            if (l_sq > r * r) t = s - q;
            else t = s + q;

            p = o + t * d;

            return IntersectionStatus.Intersect;
        }
        public static IntersectionStatus RayTriangleIntersect(Ray ray, Triangle tri, out float t, out Vector3f p)
        {
            float u, v;
            IntersectionStatus istat = RayTriangleIntersect(ray.Origin, ray.Direction, tri.V0, tri.V1, tri.V2, out u, out v, out t);
            p = tri.V[0] + u * (tri.V[1] - tri.V[0]) + v * (tri.V[2] - tri.V[0]);

            return istat;
        }

        public static IntersectionStatus RayTriangleIntersect(Ray ray, Triangle tri, out Vector3f p)
        {
            float u, v, t;
            IntersectionStatus istat = RayTriangleIntersect(ray.Origin, ray.Direction, tri.V0, tri.V1, tri.V2, out u, out v, out t);
            p = tri.V[0] + u * (tri.V[1] - tri.V[0]) + v * (tri.V[2] - tri.V[0]);

            return istat;
        }
        public static IntersectionStatus RayTriangleIntersect(Ray ray, Triangle tri, out float u, out float v, out float t)
        {
            IntersectionStatus istat = RayTriangleIntersect(ray.Origin, ray.Direction, tri.V0, tri.V1, tri.V2, out u, out v, out t);

            return istat;
        }
        public static IntersectionStatus RayTriangleIntersect(Ray ray, Triangle tri, out float u, out float v, out float t, out Vector3f p)
        {
            IntersectionStatus istat = RayTriangleIntersect(ray.Origin, ray.Direction, tri.V0, tri.V1, tri.V2, out u, out v, out t);
            p = tri.V[0] + u * (tri.V[1] - tri.V[0]) + v * (tri.V[2] - tri.V[0]);

            return istat;
        }
        public static IntersectionStatus RayTriangleIntersect(Vector3f o, Vector3f d, Vector3f v0, Vector3f v1, Vector3f v2, out float u, out float v, out float t)
        {
            // 返却する値を初期化
            u = v = t = 0.0f;

            Vector3f e1 = v1 - v0;
            Vector3f e2 = v2 - v0;

            Vector3f p = Vector3f.Cross(d, e2);
            float a = Vector3f.Dot(e1, p);

            if (a > -epsilon && a < epsilon)
            {
                u = v = t = 0.0f;
                return IntersectionStatus.Reject;
            }

            float f = 1.0f / a;
            Vector3f s = o - v0;
            u = f * Vector3f.Dot(s, p);
            if (u < 0.0 || u > 1.0)
            {
                u = v = t = 0.0f;
                return IntersectionStatus.Reject;
            }

            Vector3f q = Vector3f.Cross(s, e1);
            v = f * Vector3f.Dot(d, q);
            if (v < 0.0 || u + v > 1.0)
            {
                u = v = t = 0.0f;
                return IntersectionStatus.Reject;
            }

            t = f * Vector3f.Dot(e2, q);

            return IntersectionStatus.Intersect;
        }

        public static OverlapStatus planeBoxOverlap(Vector3f normal, Vector3f vert, Vector3f maxbox)	// -NJMP-
        {
            int q;
            Vector3f vmin = new Vector3f();
            Vector3f vmax = new Vector3f();
            float v;

            for (q = 0; q < 3; q++)
            {
                v = vert[q];					// -NJMP-

                if (normal[q] > 0.0f)
                {
                    vmin[q] = -maxbox[q] - v;	// -NJMP-
                    vmax[q] = maxbox[q] - v;	// -NJMP-
                }
                else
                {
                    vmin[q] = maxbox[q] - v;	// -NJMP-
                    vmax[q] = -maxbox[q] - v;	// -NJMP-
                }

            }

            if (Vector3f.Dot(normal, vmin) > 0.0) return  OverlapStatus.Disjoint;	// -NJMP-
            if (Vector3f.Dot(normal, vmax) >= 0.0) return  OverlapStatus.Overlap;	// -NJMP-

            return  OverlapStatus.Disjoint;

        }
        /// <summary>
        /// http://www.softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm#Line-Plane%20Intersection
        /// </summary>
        /// <param name="ray">Ray</param>
        /// <param name="pos">平面の座標</param>
        /// <param name="nor">平面の法線</param>
        /// <param name="oPos">衝突する座標</param>
        /// <returns></returns>
        public static IntersectionStatus RayPlaneIntersect(Ray ray, Vector3f pos, Vector3f nor, out Vector3f oPos)
        {
            oPos = new Vector3f();

            Vector3f u = ray.PosB - ray.PosA;
            Vector3f w = ray.PosA - pos;

            float D = Vector3f.Dot(nor, u);
            float N = -Vector3f.Dot(nor, w);

            if (Math.Abs(D) < float.Epsilon)   // segment is parallel to plane
            {
                if (N == 0)                     // segment lies in plane
                    return IntersectionStatus.Reject;
                else
                    return IntersectionStatus.Reject;   // no intersection
            }
            // they are not parallel
            // compute intersect param
            float sI = N / D;
            if (sI < 0 || sI > 1)
                return IntersectionStatus.Reject;   // no intersection

            oPos = ray.PosA + sI * u;                 // compute segment intersect point
            return IntersectionStatus.Intersect;
        }

        public static OverlapStatus TriangleAABBOverlap(AABB aabb, Triangle triangle)
        {
            Vector3f center = aabb.Center;
            Vector3f half = aabb.HalfVector;
            Vector3f[] tri = new Vector3f[3];
            tri[0] = triangle.V0;
            tri[1] = triangle.V1;
            tri[2] = triangle.V2;

            return TriangleAABBOverlap(center, half, tri);
        }

        public static OverlapStatus TriangleAABBOverlap(Vector3f boxcenter, Vector3f boxhalfsize, Vector3f[] triverts)
        {

            /*    use separating axis theorem to test overlap between triangle and box          */
            /*    need to test for overlap in these directions:                                 */
            /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle    */
            /*       we do not even need to test these)                                         */
            /*    2) normal of the triangle                                                     */
            /*    3) crossproduct(edge from tri, {x,y,z}-directin)                              */
            /*       this gives 3x3=9 more tests                                                */

            Vector3f v0, v1, v2;


            float fex, fey, fez;

            Vector3f normal, e0, e1, e2;


            /* This is the fastest branch on Sun */
            /* move everything so that the boxcenter is in (0,0,0) */
            v0 = triverts[0] - boxcenter;
            v1 = triverts[1] - boxcenter;
            v2 = triverts[2] - boxcenter;


            /* compute triangle edges */
            e0 = v1 - v0;      /* tri edge 0 */
            e1 = v2 - v1;      /* tri edge 1 */
            e2 = v0 - v2;      /* tri edge 2 */


            /* Bullet 3:  */
            /*  test the 9 tests first (this was faster) */

            fex = Math.Abs(e0.X);
            fey = Math.Abs(e0.Y);
            fez = Math.Abs(e0.Z);

            if (!AXISTEST_X01(e0.Z, e0.Y, fez, fey, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Y02(e0.Z, e0.X, fez, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Z12(e0.Y, e0.X, fey, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;



            fex = Math.Abs(e1.X);
            fey = Math.Abs(e1.Y);
            fez = Math.Abs(e1.Z);

            if (!AXISTEST_X01(e1.Z, e1.Y, fez, fey, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Y02(e1.Z, e1.X, fez, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Z0(e1.Y, e1.X, fey, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;


            fex = Math.Abs(e2.X);
            fey = Math.Abs(e2.Y);
            fez = Math.Abs(e2.Z);

            if (!AXISTEST_X2(e2.Z, e2.Y, fez, fey, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Y1(e2.Z, e2.X, fez, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;
            if (!AXISTEST_Z12(e2.Y, e2.X, fey, fex, v0, v1, v2, boxhalfsize)) return OverlapStatus.Disjoint;




            /* Bullet 1: */
            /*  first test overlap in the {x,y,z}-directions */
            /*  find min, max of the triangle each direction, and test for overlap in */
            /*  that direction -- this is equivalent to testing a minimal AABB around */
            /*  the triangle against the AABB */

            float min, max;

            /* test in X-direction */
            FINDMINMAX(v0.X, v1.X, v2.X, out min, out max);
            if (min > boxhalfsize.X || max < -boxhalfsize.X) return OverlapStatus.Disjoint;

            /* test in Y-direction */
            FINDMINMAX(v0.Y, v1.Y, v2.Y, out min, out max);
            if (min > boxhalfsize.Y || max < -boxhalfsize.Y) return OverlapStatus.Disjoint;


            /* test in Z-direction */
            FINDMINMAX(v0.Z, v1.Z, v2.Z, out min, out max);
            if (min > boxhalfsize.Z || max < -boxhalfsize.Z) return OverlapStatus.Disjoint;

            /* Bullet 2: */
            /*  test if the box intersects the plane of the triangle */
            /*  compute plane equation of triangle: normal*x+d=0 */

            normal = Vector3f.Cross(e0, e1);

            // -NJMP- (line removed here)
            if (planeBoxOverlap(normal, v0, boxhalfsize) == OverlapStatus.Disjoint) return OverlapStatus.Disjoint;	// -NJMP-

            return OverlapStatus.Overlap;   /* box and triangle overlaps */
        }

        // Real-time Rendering p.491
        public static OverlapStatus RayAABBOverlap(Ray ray, AABB aabb)
        {
            return RayAABBOverlap(ray.Center, ray.HalfVector, aabb.Center, aabb.HalfVector);
        }
        public static OverlapStatus RayAABBOverlap(Vector3f RayCenter, Vector3f RayHalf, Vector3f BoxCenter, Vector3f BoxHalf)
        {
            // BoxCenterが原点になるようにRayも平行移動
            Vector3f TransRayCenter = RayCenter - BoxCenter;

            float vx = Math.Abs(RayHalf.X);
            float vy = Math.Abs(RayHalf.Y);
            float vz = Math.Abs(RayHalf.Z);

            if (Math.Abs(TransRayCenter.X) > vx + BoxHalf.X) return OverlapStatus.Disjoint;
            if (Math.Abs(TransRayCenter.Y) > vy + BoxHalf.Y) return OverlapStatus.Disjoint;
            if (Math.Abs(TransRayCenter.Z) > vz + BoxHalf.Z) return OverlapStatus.Disjoint;

            if (Math.Abs(TransRayCenter.Y * RayHalf.Z - TransRayCenter.Z * RayHalf.Y) > BoxHalf.Y * vz + BoxHalf.Z * vy) return OverlapStatus.Disjoint;
            if (Math.Abs(TransRayCenter.X * RayHalf.Z - TransRayCenter.Z * RayHalf.X) > BoxHalf.X * vz + BoxHalf.Z * vx) return OverlapStatus.Disjoint;
            if (Math.Abs(TransRayCenter.X * RayHalf.Y - TransRayCenter.Y * RayHalf.X) > BoxHalf.X * vy + BoxHalf.Y * vx) return OverlapStatus.Disjoint;

            return OverlapStatus.Overlap;

        }

        /// <summary>
        /// [Under Construction]６枚のPlaneとRayの衝突判定（超遅い）
        /// </summary>
        /// <param name="RayCenter"></param>
        /// <param name="RayHalf"></param>
        /// <param name="BoxCenter"></param>
        /// <param name="BoxHalf"></param>
        /// <returns></returns>
        public static IntersectionStatus RayAABBIntersect_(Vector3f RayCenter, Vector3f RayHalf, Vector3f BoxCenter, Vector3f BoxHalf)
        {
            Vector3f[] p = new Vector3f[6];
            p[0] = BoxCenter - new Vector3f(BoxHalf.X, 0, 0);
            p[1] = BoxCenter + new Vector3f(BoxHalf.X, 0, 0);
            p[2] = BoxCenter - new Vector3f(0, BoxHalf.Y, 0);
            p[3] = BoxCenter + new Vector3f(0, BoxHalf.Y, 0);
            p[4] = BoxCenter - new Vector3f(0, 0, BoxHalf.Z);
            p[5] = BoxCenter + new Vector3f(0, 0, BoxHalf.Z);

            Vector3f[] n = new Vector3f[6];

            n[0] = new Vector3f(-1, 0, 0);
            n[1] = new Vector3f(1, 0, 0);
            n[2] = new Vector3f(0, -1, 0);
            n[3] = new Vector3f(0, 1, 0);
            n[4] = new Vector3f(0, 0, -1);
            n[5] = new Vector3f(0, 0, 1);

            Vector3f[] op = new Vector3f[6];
            IntersectionStatus[] istat = new IntersectionStatus[6];

            for (int i = 0; i < 6; i++)
            {
                istat[i] = RayPlaneIntersect(new Ray(RayCenter, RayHalf), p[i], n[i], out op[i]);
            }

            return IntersectionStatus.Intersect;
        }
        /// <summary>
        /// Ray-AABBの衝突判定(RTR2[Jp] p488)
        /// </summary>
        /// <param name="o"></param>
        /// <param name="d"></param>
        /// <param name="aabb"></param>
        /// <returns></returns>
        public static IntersectionStatus RayAABBIntersect(Vector3f o, Vector3f d, AABB aabb, out float t_near, out float t_far, out AABB.AABBFacet facet_near, out AABB.AABBFacet facet_far)
        {
            float t_min = float.MinValue;
            float t_max = float.MaxValue;

            t_near = 0.0f;
            t_far = 0.0f;
            facet_far = AABB.AABBFacet.NegativeX;
            facet_near = AABB.AABBFacet.NegativeX;

            Vector3f p = aabb.Center - o;


            for (int i = 0; i < 3; i++)
            {
                float e = p[i];
                float f = d[i];

                bool swap_flag = false;
                if (f * f > float.Epsilon * float.Epsilon)
                {
                    float t1 = (e + aabb.HalfVector[i]) / f;
                    float t2 = (e - aabb.HalfVector[i]) / f;
                    if (t1 > t2)
                    {
                        swap(ref t1, ref t2);
                        swap_flag = true;
                    }
                    if (t1 > t_min)
                    {
                        t_min = t1;
                        if (i == 0) {       facet_near = !swap_flag ? AABB.AABBFacet.PositiveX : AABB.AABBFacet.NegativeX; }
                        else if (i == 1) {  facet_near = !swap_flag ? AABB.AABBFacet.PositiveY : AABB.AABBFacet.NegativeY; }
                        else {              facet_near = !swap_flag ? AABB.AABBFacet.PositiveZ : AABB.AABBFacet.NegativeZ; }
                    }
                    if (t2 < t_max)
                    {
                        t_max = t2;
                        if (i == 0) {       facet_far = !swap_flag ? AABB.AABBFacet.NegativeX : AABB.AABBFacet.PositiveX; }
                        else if (i == 1) {  facet_far = !swap_flag ? AABB.AABBFacet.NegativeY : AABB.AABBFacet.PositiveY; }
                        else {              facet_far = !swap_flag ? AABB.AABBFacet.NegativeZ : AABB.AABBFacet.PositiveZ; }
                    }


                    if (t_min > t_max) return IntersectionStatus.Reject;
                    if (t_max < 0.0) return IntersectionStatus.Reject;
                }
                else if (-e - aabb.HalfVector[i] > 0.0 || -e + aabb.HalfVector[i] < 0.0) return IntersectionStatus.Reject;
            }


            if (t_min > 0.0) 
            {
                // AABBの外側に始点があるとき
                t_near = t_min;
                t_far = t_max;
                return IntersectionStatus.Intersect;
            }
            else 
            {
                // AABBの内側に始点があるとき
                t_near = t_min;
                t_far = t_max;
                return IntersectionStatus.Intersect;
            }

        }
        private static void swap(ref float x, ref float y)
        {
            float t = x;
            x = y;
            y = t;
        }
        public static OverlapStatus AABBOverlap(AABB a, AABB b)
        {
            Vector3f MinA = a.Center - a.HalfVector;
            Vector3f MaxA = a.Center + a.HalfVector;
            Vector3f MinB = b.Center - b.HalfVector;
            Vector3f MaxB = b.Center + b.HalfVector;
            for (int i = 0; i < 3; i++)
            {
                if (MinA[i] > MaxB[i] || MinB[i] > MaxA[i]) return OverlapStatus.Disjoint;
            }

            return OverlapStatus.Overlap;
        }
        public static OverlapStatus PointAABBOverlap(Vector3f p, AABB aabb)
        {
            Vector3f q = p - aabb.Center;

            if (q.X < -aabb.HalfVector.X || aabb.HalfVector.X < q.X) return OverlapStatus.Disjoint;
            if (q.Y < -aabb.HalfVector.Y || aabb.HalfVector.Y < q.Y) return OverlapStatus.Disjoint;
            if (q.Z < -aabb.HalfVector.Z || aabb.HalfVector.Z < q.Z) return OverlapStatus.Disjoint;

            return OverlapStatus.Overlap;
        }
        // デバッグしてない
        // リアルタイムレンダリング 第２版（日本語版） pp.509
        public static OverlapStatus SphereAABBOverlap(Vector3f center, float radius, AABB aabb)
        {
            float d = 0;
            Vector3f aabb_min = aabb.Center - aabb.HalfVector;
            Vector3f aabb_max = aabb.Center + aabb.HalfVector;

            if      (center.X < aabb_min.X) d = d + (center.X - aabb_min.X) * (center.X - aabb_min.X);
            else if (center.X > aabb_max.X) d = d + (center.X - aabb_max.X) * (center.X - aabb_max.X);
            if      (center.Y < aabb_min.Y) d = d + (center.Y - aabb_min.Y) * (center.Y - aabb_min.Y);
            else if (center.Y > aabb_max.Y) d = d + (center.Y - aabb_max.Y) * (center.Y - aabb_max.Y);
            if      (center.Z < aabb_min.Z) d = d + (center.Z - aabb_min.Z) * (center.Z - aabb_min.Z);
            else if (center.Z > aabb_max.Z) d = d + (center.Z - aabb_max.Z) * (center.Z - aabb_max.Z);

            if (d > radius * radius) return OverlapStatus.Disjoint;
            else return OverlapStatus.Overlap;
        }
        //static int triBoxOverlap(Vector3d boxcenter, Vector3d boxhalfsize, Triangle triangle)
        //{



        //    /*    use separating axis theorem to test overlap between triangle and box */

        //    /*    need to test for overlap in these directions: */

        //    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */

        //    /*       we do not even need to test these) */

        //    /*    2) normal of the triangle */

        //    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */

        //    /*       this gives 3x3=9 more tests */

        //    Vector3d v0, v1, v2;

        //    //   float axis[3];

        //    float min, max, p0, p1, p2, rad, fex, fey, fez;		// -NJMP- "d" local variable removed

        //    Vector3d normal, e0, e1, e2;



        //    /* This is the fastest branch on Sun */

        //    /* move everything so that the boxcenter is in (0,0,0) */

        //    v0 = triangle.V0 - boxcenter;

        //    v1 = triangle.V1 - boxcenter;

        //    v2 = triangle.V2 - boxcenter;



        //    /* compute triangle edges */

        //    e0 = v1 - v0;      /* tri edge 0 */

        //    e1 = v2 - v1;      /* tri edge 1 */

        //    e2 = v0 - v2;      /* tri edge 2 */



        //    /* Bullet 3:  */

        //    /*  test the 9 tests first (this was faster) */

        //    fex = Math.Abs(e0.X);

        //    fey = Math.Abs(e0.Y);

        //    fez = Math.Abs(e0.Z);

        //    AXISTEST_X01(e0.Z, e0.Y, fez, fey, v0, v1, v2, boxhalfsize);
        //    AXISTEST_Y02(e0.Z, e0.X, fez, fex, v0, v1, v2, boxhalfsize);
        //    AXISTEST_Z12(e0.Y, e0.X, fey, fex, v0, v1, v2, boxhalfsize);

        //    fex = Math.Abs(e1.X);

        //    fey = Math.Abs(e1.Y);

        //    fez = Math.Abs(e1.Z);

        //    AXISTEST_X01(e1.Z, e1.Y, fez, fey, v0, v1, v2, boxhalfsize);
        //    AXISTEST_Y02(e1.Z, e1.X, fez, fex, v0, v1, v2, boxhalfsize);
        //    AXISTEST_Z0(e1.Y, e1.X, fey, fex, v0, v1, v2, boxhalfsize);



        //    fex = Math.Abs(e2.X);
        //    fey = Math.Abs(e2.Y);
        //    fez = Math.Abs(e2.Z);

        //    AXISTEST_X2(e2.Z, e2.Y, fez, fey, v0, v1, v2, boxhalfsize);


        //    AXISTEST_Y1(e2.Z, e2.X, fez, fex, v0, v1, v2, boxhalfsize);


        //    AXISTEST_Z12(e2.Y, e2.X, fey, fex, v0, v1, v2, boxhalfsize);




        //    /* Bullet 1: */

        //    /*  first test overlap in the {x,y,z}-directions */

        //    /*  find min, max of the triangle each direction, and test for overlap in */

        //    /*  that direction -- this is equivalent to testing a minimal AABB around */

        //    /*  the triangle against the AABB */



        //    /* test in X-direction */

        //    FINDMINMAX(v0.X, v1.X, v2.X, out min, out max);


        //    if (min > boxhalfsize.X || max < -boxhalfsize.X) return 0;



        //    /* test in Y-direction */

        //    FINDMINMAX(v0.Y, v1.Y, v2.Y, out min, out max);

        //    if (min > boxhalfsize.Y || max < -boxhalfsize.Y) return 0;



        //    /* test in Z-direction */

        //    FINDMINMAX(v0.Z, v1.Z, v2.Z, out min, out max);

        //    if (min > boxhalfsize.Z || max < -boxhalfsize.Z) return 0;



        //    /* Bullet 2: */

        //    /*  test if the box intersects the plane of the triangle */

        //    /*  compute plane equation of triangle: normal*x+d=0 */

        //    normal = Vector3d.Cross(e0, e1);

        //    // -NJMP- (line removed here)

        //    //if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;	// -NJMP-

        //    //if(planeBoxOverlap(normal, v0, boxhalfsize) != 0)



        //    return 1;   /* box and triangle overlaps */

        //}
    }
}
