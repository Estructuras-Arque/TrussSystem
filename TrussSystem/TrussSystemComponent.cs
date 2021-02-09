using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;
using Rhino.Geometry.Intersect;

// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace TrussSystem
{
    public class TrussSystemComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public TrussSystemComponent()
          : base("Portic/Truss System", "TPS",
              "Description",
              "Arque Structures", "Portics")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPlaneParameter("plane", "pl", "description", GH_ParamAccess.item, Plane.WorldXY);
            pManager.AddIntegerParameter("typology", "ty", "description", GH_ParamAccess.item, 3);
            pManager.AddIntegerParameter("base type", "bT", "description", GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("support type", "sT", "description", GH_ParamAccess.item, 1);
            pManager.AddIntegerParameter("span left", "sL", "description", GH_ParamAccess.item, 5000);
            pManager.AddIntegerParameter("span right", "sR", "description", GH_ParamAccess.item, 5000);
            pManager.AddIntegerParameter("max height", "mH", "description", GH_ParamAccess.item, 3500);
            pManager.AddIntegerParameter("clear height", "cH", "description", GH_ParamAccess.item, 2700);
            pManager.AddIntegerParameter("right height", "rH", "description", GH_ParamAccess.item, 3000);
            pManager.AddIntegerParameter("left height", "lH", "description", GH_ParamAccess.item, 3000);
            pManager.AddIntegerParameter("subdivision count", "sC", "description", GH_ParamAccess.item, 4);
            pManager.AddIntegerParameter("min length", "miL", "description", GH_ParamAccess.item, 1500);
            pManager.AddIntegerParameter("max length", "maL", "description", GH_ParamAccess.item, 2000);
            pManager.AddIntegerParameter("portic type", "pT", "description", GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("truss type", "tT", "description", GH_ParamAccess.item, 0);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("columns", "cl", "description", GH_ParamAccess.list);
            pManager.AddCurveParameter("superior curves", "sC", "description", GH_ParamAccess.list);
            pManager.AddPointParameter("superior points", "sP", "description", GH_ParamAccess.list);
            pManager.AddCurveParameter("inferior curves", "iC", "description", GH_ParamAccess.list);
            pManager.AddPointParameter("inferior points", "iP", "description", GH_ParamAccess.list);
            pManager.AddCurveParameter("intermediate curves", "iC", "description", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int porticType = 0;
            int trussType = 0;
            int typology = 0;
            int clHeight = 0;
            int crHeight = 0;
            int maxHeight = 0;
            int spanOne = 0;
            int spanTwo = 0;
            int clearHeight = 0;
            int baseType = 0;
            int subdCount = 0;
            int supType = 0;
            int minLength = 0;
            int maxLength = 0;
            Plane inputPlane = new Plane();
            if (!DA.GetData(0, ref inputPlane)) return;
            if (!DA.GetData(1, ref typology)) return;
            if (!DA.GetData(2, ref baseType)) return;
            if (!DA.GetData(3, ref supType)) return;
            if (!DA.GetData(4, ref spanOne)) return;
            if (!DA.GetData(5, ref spanTwo)) return;
            if (!DA.GetData(6, ref maxHeight)) return;
            if (!DA.GetData(7, ref clearHeight)) return;
            if (!DA.GetData(8, ref crHeight)) return;
            if (!DA.GetData(9, ref clHeight)) return;
            if (!DA.GetData(10, ref subdCount)) return;
            if (!DA.GetData(11, ref minLength)) return;
            if (!DA.GetData(12, ref maxLength)) return;
            if (!DA.GetData(13, ref porticType)) return;
            if (!DA.GetData(14, ref trussType)) return;

            // initialize outputs
            List<Point3d> superiorBasePoints = new List<Point3d>();
            List<Point3d> inferiorBasePoints = new List<Point3d>();


            List<Curve> columnCurves = new List<Curve>();
            List<Curve> superiorBaseCurves = new List<Curve>();
            List<Curve> inferiorBaseCurves = new List<Curve>();
            List<Curve> intermediateCurvesList = new List<Curve>();

            DataTree<Point3d> superiorPoints = new DataTree<Point3d>();
            DataTree<Point3d> inferiorPoints = new DataTree<Point3d>();
            DataTree<Curve> superiorCurves = new DataTree<Curve>();
            DataTree<Curve> inferiorCurves = new DataTree<Curve>();
            DataTree<Curve> intermediateCurves= new DataTree<Curve>();
            DataTree<Point3d> thickPoints = new DataTree<Point3d>();
            DataTree<Point3d> straightPoints = new DataTree<Point3d>();


            List<Vector3d> sideVectors = new List<Vector3d>();
            List<double> hypotenusVariables = new List<double>();

            // warren division restrict
            if (trussType == 2 && subdCount % 2 == 1) subdCount += 1;

            // get minimum height and its index to select the rest of the elements according to it
            int[] columnHeights = new int[] { clHeight, crHeight };
             // Print("height left: {0} height right: {1}", columnHeights[0], columnHeights[1]);
             int cMinHeight = columnHeights.Min();
             int cMaxHeight = columnHeights.Max();

            // restrict min height to max column height + 200 if typology is curved or piched
            if(maxHeight <= cMaxHeight && typology == 1) { maxHeight = cMaxHeight + 200; }
            else if (maxHeight <= cMaxHeight && (typology == 2 || typology == 3))
            {
               maxHeight = cMaxHeight;
            }

            // get min and max heights indices
             int notSelected = Array.IndexOf(columnHeights, cMaxHeight);
            int index = Array.IndexOf(columnHeights, cMinHeight);

            // generate superior truss base points
            superiorBasePoints = SuperiorBasePoints(inputPlane,typology, spanOne, spanTwo, maxHeight, ref crHeight, ref clHeight, cMinHeight, cMaxHeight);

            // generate the columns
            columnCurves = ColumnCurves(inputPlane, typology, spanOne, spanTwo, superiorBasePoints, columnCurves);

            // generate superior truss curve
            superiorBaseCurves = GenerateThickCurves(typology, maxHeight, crHeight, clHeight, superiorBasePoints, 0);

            // calculate the clear difference , and if the clear height exceeds minimum column height, restrict it to 200
            // Print("smallest height: {0} index: {1}", cMinHeight, index);
            int difference = ComputeClearHeight(typology, maxHeight, ref clearHeight, cMinHeight);

            // ====== generate inferior truss points & curves ======= //

            // compute the offset variable
            double offsetFactor = ComputeOffset(index, difference, superiorBaseCurves[index]);

            // compute and store the normal vectors at each point
            List<Vector3d> normalVectors = GetFirstLastNormalVectors(typology, superiorBaseCurves);
            sideVectors.Add(normalVectors[0]);
            sideVectors.Add(normalVectors[2]);
            normalVectors[1].Rotate(-0.5 * Math.PI, Vector3d.YAxis);

            // compute highest hypotenus and middle hypotenus
            double middleHypotenus = ComputeCenterHypotenus(offsetFactor, superiorBaseCurves[index], normalVectors[1]);
            double hypotenus2 = ComputeCornerHypotenus(notSelected, offsetFactor, superiorBaseCurves[notSelected]);


            // store points at normal vectors at each point
            for (int i = 0; i < normalVectors.Count; i++)
            {
                Point3d tempPt = new Point3d(i != 1 ? offsetFactor * normalVectors[i] + superiorBasePoints[i] : middleHypotenus * normalVectors[i] + superiorBasePoints[i]);
                inferiorBasePoints.Add(tempPt);
            }

            // on each side create the clear point at its corresponding clear height
            for (int i = 0; i < sideVectors.Count; i++)
            {
                Point3d tempPt = new Point3d(i != index ? -hypotenus2 * Vector3d.ZAxis + superiorBasePoints[i == 1 ? 2 : 0] : -difference * Vector3d.ZAxis + superiorBasePoints[i == 1 ? 2 : 0]);
                inferiorBasePoints.Insert(i == 0 ? 0 : 4, tempPt);
            }
            List<Point3d> rigidPoints = new List<Point3d>();
            Point3d clearPoint = inferiorBasePoints[index == 0 ? 0 : 4];
            rigidPoints.Add(inferiorBasePoints[0]);
            Point3d counterClearPoint = inferiorBasePoints[notSelected == 0 ? 0 : 4];
            rigidPoints.Add(inferiorBasePoints[4]);
            // Print("count {0}", inferiorBasePoints.Count);
            // generate lower curves
            List<Curve> thickCurvesList = GenerateThickCurves(typology, maxHeight, crHeight, clHeight, inferiorBasePoints, 1);

            List<Curve> straightCurvesList = GenerateStraightCurves(clearPoint, counterClearPoint, superiorBasePoints);

            if (baseType == 1) inferiorBaseCurves = thickCurvesList;
            else inferiorBaseCurves = straightCurvesList;
            // Print("count {0}", inferiorBaseCurves.Count);

            // compute subdivision count

            List<int> computedDivisions = ComputeDivisions(superiorBaseCurves, subdCount, trussType, minLength, maxLength);


            // iterate on each base curve and divide by count
            for (int i = 0; i < superiorBaseCurves.Count; i++)
            {
                Point3d[] upperSubdPoints = new Point3d[] { };
                List<Point3d> tempIntersections = new List<Point3d>();

                List<Point3d> thickPointsList = new List<Point3d>();
                List<Point3d> straightPointsList = new List<Point3d>();
                Curve top_crv = superiorBaseCurves[i];
                Curve bottom_crv = inferiorBaseCurves[i];

                GH_Path path = new GH_Path(i);
                // store parameters t on top_crv for each point
                double[] subdParameters = top_crv.DivideByCount(computedDivisions[i], true, out upperSubdPoints);
                for (int j = 0; j < upperSubdPoints.Length; j++)
                {
                    if (trussType <= 1 || baseType == 0)
                    {
                        Plane yzPlane = new Plane(upperSubdPoints[j], Vector3d.ZAxis, Vector3d.YAxis);
                        var events = Intersection.CurvePlane(Curve.JoinCurves(inferiorBaseCurves)[0], yzPlane, 0.001);
                        if (events != null)
                        {
                            IntersectionEvent cpx_event = new IntersectionEvent();
                            for (int k = 0; k < events.Count; k++)
                            {
                                cpx_event = events[k];
                            }
                            if (j == upperSubdPoints.Length - 1 && typology == 3 && baseType == 1)
                            {
                                tempIntersections.Add(inferiorBasePoints[2]);
                            }
                            else
                            {
                                tempIntersections.Add(cpx_event.PointA);
                            }
                        }
                    }
                    else
                    {
                        double t;
                        bottom_crv.ClosestPoint(upperSubdPoints[j], out t, 0);
                        Point3d tempPoint = bottom_crv.PointAt(t);
                        Point3d endPt = bottom_crv.PointAtLength(bottom_crv.GetLength());
                        if (j != upperSubdPoints.Length - 1) tempIntersections.Add(tempPoint);
                        else tempIntersections.Add(endPt);

                    }
                    if (baseType == 1) thickPointsList = tempIntersections;
                    else straightPointsList = tempIntersections;
                }
                superiorPoints.AddRange(upperSubdPoints, path);
                thickPoints.AddRange(thickPointsList, path);
                straightPoints.AddRange(straightPointsList, path);
            }

            // generate straight inferior Points
            if (baseType == 0)
            {
                inferiorPoints = straightPoints;
                inferiorBaseCurves = straightCurvesList;
            }
            // generate thickened inferior points
            else
            {
                inferiorPoints = thickPoints;
                inferiorBaseCurves = thickCurvesList;
            }
            List<Point3d> tempSuperiorPoints = new List<Point3d>();
            List<Point3d> tempInferiorPoints = new List<Point3d>();
            List<Curve> superiorCurvesList = new List<Curve>();
            List<Curve> inferiorCurvesList = new List<Curve>();

            // get intermediate connections by truss type
            for (int i = 0; i < superiorBaseCurves.Count; i++)
            {
                tempSuperiorPoints = superiorPoints.Branch(i);
                tempInferiorPoints = inferiorPoints.Branch(i);
                List<Curve> tempIntermediateList = new List<Curve>();

                for (int j = 0; j < tempSuperiorPoints.Count; j++)
                {
                    int[] path = new int[] { i, j };
                    Line line1 = new Line();
                    Line line2 = new Line();
                    Line line3 = new Line();
                    // HOWE
                    if (trussType == 0)
                    {
                        if (supType == 1)
                        {
                            if (j != 0) line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j - 1]);
                            if (j != 0) line2 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j]);
                        }
                        else
                        {
                            if (j > 1) line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j - 1]);
                            if (j == 0) line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j + 1]);
                            if (j != 0) line2 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j]);
                        }

                    }

                    // PRATT
                    else if (trussType == 1)
                    {
                        if (j != 0) line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j]);
                        if (j != tempSuperiorPoints.Count - 1)
                        {
                            line2 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j + 1]);
                        }
                    }

                    // WARREN
                    else if (trussType == 2 && j % 2 == 0)
                    {
                        if (j != tempSuperiorPoints.Count - 1) line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j + 1]);
                        if (j != 0) line2 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j - 1]);
                    }

                    // WARREN INTERLOCK
                    else if (trussType == 3)
                    {
                        if (j % 2 == 0)
                        {
                            if (j != tempSuperiorPoints.Count - 1)
                            {
                                line1 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j + 1]);
                            }
                            if (j != 0) line3 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j - 1]);
                        }
                        else
                        {
                            line2 = new Line(tempSuperiorPoints[j], tempInferiorPoints[j]);
                        }
                    }
                    if (line1.Length != 0) tempIntermediateList.Add(line1.ToNurbsCurve());
                    if (line2.Length != 0) tempIntermediateList.Add(line2.ToNurbsCurve());
                    if (trussType == 3 && line3.Length != 0) tempIntermediateList.Add(line3.ToNurbsCurve());
                }
                intermediateCurves.AddRange(tempIntermediateList);
            }
            intermediateCurvesList = intermediateCurves.Branch(0);
            //DataTree<Point3d> quotationPoints = QuotationPoints(superiorBasePoints);

            // get superior and inferior connections by truss type

            superiorPoints.Branch(1).Reverse();
            superiorPoints.Branch(1).RemoveAt(0);
            superiorPoints.Flatten(null);
            List<Point3d> superiorPointsList = superiorPoints.Branch(0);
            // Print("{0}", superiorPointsList.Count);
            inferiorPoints.Branch(1).Reverse();
            inferiorPoints.Branch(1).RemoveAt(0);
            inferiorPoints.Flatten(null);
            List<Point3d> inferiorPointsList = inferiorPoints.Branch(0);
            // Print("{0}", inferiorPointsList.Count);
            if (trussType > 1 && baseType == 1)
            {
                inferiorPointsList.Insert(0, inferiorBasePoints[0]);
                inferiorPointsList.Insert(inferiorPointsList.Count, inferiorBasePoints[4]);
            }

            //for (int i = 0; i < superiorPointsList.Count; i++)
            //{
            //    Line supCurve = new Line();

            //    if (trussType != 2)
            //    {
            //        if (i != superiorPointsList.Count - 1)
            //        {
            //            supCurve = new Line(superiorPointsList[i], superiorPointsList[i + 1]);
            //        }
            //    }
            //    else if (trussType == 2)
            //    {

            //        if (i % 2 == 0 && i != superiorPointsList.Count - 1)
            //        {
            //            supCurve = new Line(superiorPointsList[i], superiorPointsList[i + 2]);
            //        }

            //    }
            //    if (supCurve.IsValid && supCurve.Length != 0)
            //    {
            //        superiorCurvesList.Add(supCurve.ToNurbsCurve());
            //    }
            //}

            //for (int i = 0; i < inferiorPointsList.Count; i++)
            //{
            //    Line infCurve = new Line();

            //    if (trussType <= 1)
            //    {
            //        if (i != inferiorPointsList.Count - 1)
            //        {
            //            infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 1]);
            //        }
            //    }
            //    else if (trussType == 2)
            //    {
            //        if (baseType == 0)
            //        {
            //            if (i < 1) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 1]);
            //            else if (i == inferiorPointsList.Count - 2) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 1]);
            //            else if (i % 2 == 1)
            //            {
            //                infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 2]);
            //            }
            //        }
            //        else
            //        {
            //            if (i % 2 == 0)
            //            {
            //                if (i != inferiorPointsList.Count - 1) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 2]);
            //            }
            //        }
            //    }
            //    else if (trussType == 3)
            //    {
            //        if (baseType == 0)
            //        {
            //            if (i % 2 == 1 && i < inferiorPointsList.Count - 2) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 2]);
            //            else if (i <= 0 || i == inferiorPointsList.Count - 2) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 1]);
            //        }
            //        else
            //        {
            //            if (i % 2 == 0 && i < inferiorPointsList.Count - 1) infCurve = new Line(inferiorPointsList[i], inferiorPointsList[i + 2]);
            //        }
            //    }
            //    if (infCurve.IsValid && infCurve.Length != 0)
            //    {
            //        inferiorCurvesList.Add(infCurve.ToNurbsCurve());
            //    }
            //}
            if (supType == 0)
            {
                List<Curve> tempList = inferiorBaseCurves;
                //inferiorBaseCurves = new List<Curve>();
                //inferiorCurvesList.RemoveAt(0);
                //inferiorCurvesList.RemoveAt(inferiorCurvesList.Count - 1);
                //for (int i = 0; i < tempList.Count; i++)
                //{
                //    double t;
                //    tempList[i].ClosestPoint(inferiorPointsList[i == 0 ? 1 : inferiorPointsList.Count - 2], out t);
                //    inferiorBaseCurves.Add(tempList[i].Split(t)[1]);
                //}
                //inferiorPointsList.RemoveAt(0);
                //inferiorPointsList.RemoveAt(inferiorPointsList.Count - 1);

            }

            if (porticType == 1)
            {
                DA.SetDataList(0, columnCurves);
                DA.SetDataList(1, superiorBaseCurves);
                DA.SetDataList(2, superiorPointsList);
                DA.SetDataList(3, inferiorBaseCurves);
                DA.SetDataList(4, inferiorPointsList);
                DA.SetDataList(5, intermediateCurvesList);
            }
            else
            {
                DA.SetDataList(0, columnCurves);
                DA.SetDataList(1, superiorBaseCurves);
                DA.SetDataList(2, superiorPointsList);

            }
        }

        private List<int> ComputeDivisions(List<Curve> superiorBaseCurves, int subdCount, int trussType, int minLength, int maxLength)
        {
            List<int> divisions = new List<int>();
            double panelLength;
            double barLength;

            Interval lengthsInterval = new Interval(minLength, maxLength);
            for (int i = 0; i < superiorBaseCurves.Count; i++)
            {
                int count;
                barLength = superiorBaseCurves[i].GetLength();
                panelLength = barLength / subdCount;
                int minDivision = Convert.ToInt32(barLength / lengthsInterval.T0);

                int maxDivision = Convert.ToInt32(barLength / lengthsInterval.T1);
                if (lengthsInterval.IncludesParameter(panelLength))
                {
                    count = subdCount;
                    if (trussType == 2) count *= 2;
                    divisions.Add(count);
                }
                else if (panelLength < lengthsInterval.T0)
                {
                    count =Convert.ToInt32( minDivision);
                    if (trussType == 2) count *= 2;

                    divisions.Add(count);
                }
                else if (panelLength > lengthsInterval.T1)
                {
                    count = Convert.ToInt32(maxDivision);
                    if (trussType == 2) count *= 2;

                    divisions.Add(count);
                }
            }
            for (int i = 0; i < divisions.Count; i++)
            {
                if(trussType >= 2 && divisions[i]%2==1)
                {
                    int tempDiv = divisions[i] + 1;
                    divisions.RemoveAt(i);
                    divisions.Insert(i, tempDiv);
                }

            }
            return divisions;
        }

        private List<Curve> GenerateStraightCurves(Point3d clearPoint, Point3d counterClearPoint, List<Point3d> superiorBasePoints)
        {
            List<Curve> straightCurvesList = new List<Curve>();
            Point3d pt1 = new Point3d(superiorBasePoints[0].X, superiorBasePoints[0].Y, clearPoint.Z);
            Point3d middlePoint = new Point3d(superiorBasePoints[1].X, superiorBasePoints[1].Y, clearPoint.Z);
            Point3d pt2 = new Point3d(superiorBasePoints[2].X, superiorBasePoints[2].Y, clearPoint.Z);
            Line line1 = new Line(pt1, middlePoint);
            straightCurvesList.Add(line1.ToNurbsCurve());
            Line line2 = new Line(pt2, middlePoint);
            straightCurvesList.Add(line2.ToNurbsCurve());
            return straightCurvesList;
        }

        private double ComputeCornerHypotenus(int index, double offsetFactor, Curve highestCurve)
        {
            // get tangent vector at start point
            Vector3d tanVecStart = highestCurve.TangentAtStart;

            // compute at smallest index, the move factor => sin(angle) * hypotenus = height
            double angle = new double();
            if (index == 0) angle = Vector3d.VectorAngle(tanVecStart, Vector3d.XAxis);
            else angle = Vector3d.VectorAngle(tanVecStart, -Vector3d.XAxis);
            double theta = 0.5 * (Math.PI) - angle;
            double moveFactor = offsetFactor / Math.Sin(theta);
            return moveFactor;
        }

        private double ComputeCenterHypotenus(double offsetFactor, Curve curve, Vector3d normalVector)
        {
            // get tangent vector at start point
            Vector3d tanVecEnd = curve.TangentAtEnd;
            // compute at smallest index, the move factor => sin(angle) * hypotenus = height
            double angle = Vector3d.VectorAngle(tanVecEnd, normalVector);
            double moveFactor = offsetFactor / Math.Sin(angle);
            return moveFactor;
        }

        private List<Vector3d> GetFirstLastNormalVectors(int typology, List<Curve> superiorBaseCurves)
        {
            List<Vector3d> tangentVectors = new List<Vector3d>();
            Curve[] joinCurves = Curve.JoinCurves(superiorBaseCurves);
            Curve superiorCurve = joinCurves[0];
            Vector3d tangentAtStart = new Vector3d(typology == 1 ? superiorCurve.TangentAtStart : superiorBaseCurves[0].PointAtEnd - superiorBaseCurves[0].PointAtStart);
            tangentAtStart.Unitize();
            Vector3d tangentAtEnd = new Vector3d(typology == 1 ? superiorCurve.TangentAtEnd : superiorBaseCurves[1].PointAtEnd - superiorBaseCurves[1].PointAtStart);
            tangentAtEnd.Unitize();
            Vector3d tangent1 = new Vector3d(superiorBaseCurves[0].PointAtStart - superiorBaseCurves[0].PointAtEnd);
            tangent1.Unitize();
            Vector3d tangent2 = new Vector3d(superiorBaseCurves[1].PointAtStart - superiorBaseCurves[1].PointAtEnd);
            tangent2.Unitize();
            double centerAngle = Vector3d.VectorAngle(tangent1, tangent2);
            tangent1.Rotate(-centerAngle / 2, Vector3d.YAxis);
            tangent2.Rotate(centerAngle / 2, Vector3d.YAxis);
            Vector3d tangentCenter = new Vector3d(tangent1 + tangent2 / 2);
            tangentCenter.Unitize();
            tangentAtStart.Rotate(0.5 * Math.PI, Vector3d.YAxis);
            tangentAtEnd.Rotate(typology == 1 ? 0.5 * Math.PI : -0.5 * Math.PI, Vector3d.YAxis);
            tangentVectors.Add(tangentAtStart);
            tangentVectors.Add(tangentCenter);
            tangentVectors.Add(tangentAtEnd);
            return tangentVectors;

        }

        private double ComputeOffset(int index, int difference, Curve lowestCurve)
        {
            // get tangent vector at start point
            Vector3d tanVecStart = lowestCurve.TangentAtStart;

            // compute at smallest index, the move factor => sin(angle) * hypotenus = height
            double angle = new double();
            if (index == 0) angle = Vector3d.VectorAngle(tanVecStart, Vector3d.XAxis);
            else angle = Vector3d.VectorAngle(tanVecStart, -Vector3d.XAxis);
            double theta = 0.5 * (Math.PI) - angle;
            double moveFactor = Math.Sin(theta) * difference;
            return moveFactor;
        }

        private int ComputeClearHeight(int typology, int maxHeight, ref int clearHeight, int cMinHeight)
        {
            if (typology == 0) cMinHeight = maxHeight;
            int difference = new int();
            if (clearHeight >= cMinHeight - 200)
            {
                difference = 200;
            }
            else
            {
                difference = cMinHeight - clearHeight;
            }
            return difference;
        }

        private List<Curve> GenerateThickCurves(int typology, int maxHeight, int crHeight, int clHeight, List<Point3d> basePointsList, int level)
        {
            List<Curve> baseCurves = new List<Curve>();
            Vector3d vector1 = new Vector3d(basePointsList[1] - basePointsList[0]);
            Vector3d vector2 = new Vector3d(basePointsList[1] - basePointsList[2]);
            double angle = Vector3d.VectorAngle(vector1, vector2);

            // draw straight lines
            if ((typology != 1) || (maxHeight == crHeight && maxHeight == clHeight || angle == 0 || angle == Math.PI))
            {
                // draw the typology curve
                Line line1 = new Line(basePointsList[0], basePointsList[level == 0 ? 1 : 2]);
                Line line2 = new Line(basePointsList[level == 0 ? 2 : 4], basePointsList[level == 0 ? 1 : 2]);
                baseCurves.Add(line1.ToNurbsCurve());
                baseCurves.Add(line2.ToNurbsCurve());
            }

            // else draw arch
            else
            {
                Arc arc = new Arc(basePointsList[level == 0 ? 0 : 1], basePointsList[level == 0 ? 1 : 2], basePointsList[level == 0 ? 2 : 3]);
                Point3d pt1 = new Point3d((basePointsList[level == 0 ? 0 : 1].X + basePointsList[level == 0 ? 1 : 2].X) / 2, (basePointsList[level == 0 ? 0 : 1].Y + basePointsList[level == 0 ? 1 : 2].Y) / 2, (basePointsList[level == 0 ? 0 : 1].Z + basePointsList[level == 0 ? 1 : 2].Z) / 2);
                Point3d pt2 = new Point3d((basePointsList[level == 0 ? 1 : 2].X + basePointsList[level == 0 ? 2 : 3].X) / 2, (basePointsList[level == 0 ? 1 : 2].Y + basePointsList[level == 0 ? 2 : 3].Y) / 2, (basePointsList[level == 0 ? 1 : 2].Z + basePointsList[level == 0 ? 2 : 3].Z) / 2);
                pt1 = arc.ClosestPoint(pt1);
                pt2 = arc.ClosestPoint(pt2);
                Arc arc1 = new Arc(basePointsList[level == 0 ? 0 : 1], pt1, basePointsList[level == 0 ? 1 : 2]);
                Arc arc2 = new Arc(basePointsList[level == 0 ? 2 : 3], pt2, basePointsList[level == 0 ? 1 : 2]);

                if (level == 1)
                {
                    List<Curve> tempList = new List<Curve>();
                    Line tangentStart = new Line(basePointsList[0], basePointsList[1]);
                    tempList.Add(tangentStart.ToNurbsCurve());
                    tempList.Add(arc1.ToNurbsCurve());
                    baseCurves.Add(Curve.JoinCurves(tempList, 0.001, true)[0]);
                    tempList = new List<Curve>();
                    Line tangentEnd = new Line(basePointsList[4], basePointsList[3]);
                    tempList.Add(tangentEnd.ToNurbsCurve());
                    tempList.Add(arc2.ToNurbsCurve());
                    baseCurves.Add(Curve.JoinCurves(tempList, 0.001, true)[0]);
                }
                else
                {
                    baseCurves.Add(arc1.ToNurbsCurve());
                    baseCurves.Add(arc2.ToNurbsCurve());
                }

            }

            return baseCurves;
        }

        private List<Curve> ColumnCurves(Plane plane,int typology, int spanOne, int spanTwo, List<Point3d> superiorBasePoints, List<Curve> columnCurves)
        {
            Point3d colPt1 = new Point3d(typology == 3 ? plane.Origin.X-spanOne : plane.Origin.X-spanOne / 2, plane.Origin.Y, 0);
            Point3d colPt2 = new Point3d(typology == 3 ? plane.Origin.X+spanTwo : plane.Origin.X+spanOne / 2, plane.Origin.Y, 0);
            columnCurves.Add(new Line(colPt1, superiorBasePoints[0]).ToNurbsCurve());
            columnCurves.Add(new Line(colPt2, superiorBasePoints[2]).ToNurbsCurve());
            return columnCurves;
        }

        private List<Point3d> SuperiorBasePoints(Plane plane, int typology, int spanOne, int spanTwo, int maxHeight, ref int crHeight, ref int clHeight, int cMinHeight, int cMaxHeight)
        {
            
            if (typology == 1)
            {
                if (maxHeight > spanOne * 0.8) maxHeight = Convert.ToInt32(spanOne * 0.8);
                if (crHeight > spanOne * 0.7) crHeight = Convert.ToInt32(spanOne * 0.7);
                if (clHeight > spanOne * 0.7) clHeight = Convert.ToInt32(spanOne * 0.7);
            }
            List<Point3d> superiorBasePoints = new List<Point3d>();
            Point3d pt1 = new Point3d(typology == 3 ? plane.Origin.X-spanOne : plane.Origin.X - spanOne / 2, plane.Origin.Y, typology == 0 ? maxHeight : clHeight);
            Point3d pt3 = new Point3d(typology == 3 ? plane.Origin.X+spanTwo : plane.Origin.X+spanOne / 2, plane.Origin.Y, typology == 0 ? maxHeight : crHeight);
            Point3d pt2 = new Point3d((pt3.X + pt1.X) / 2, (pt3.Y + pt1.Y) / 2, typology == 2 ? (pt3.Z + pt1.Z) / 2 : maxHeight);
            superiorBasePoints.Add(pt1);
            superiorBasePoints.Add(pt2);
            superiorBasePoints.Add(pt3);
            return superiorBasePoints;
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return Properties.Resources.portic;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("15702a80-a0b6-4df2-9db3-aacb1deecef6"); }
        }
    }
}
