#include "vtkChartXYZ.h"
#include "vtkContextMouseEvent.h"
#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPlotLine3D.h"
#include "vtkPlotPoints3D.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkTable.h"
#include "vtkUnsignedCharArray.h"
#include "vtkVector.h"
#include "vtkPen.h"
#include "vtkMath.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkPlotPoints3D.h"
#include "vtkPlotLine3D.h"
#include "vtkDataSetAttributes.h"
#include <vtkOutputWindow.h>
#include <vtkSmartPointer.h>

#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>


struct Point3D {
    double x, y, z;
};

Point3D generateRandomPointAround(const Point3D& center, const Point3D& variance) {
    Point3D p;
    p.x = center.x + vtkMath::Gaussian(0, variance.x);
    p.y = center.y + vtkMath::Gaussian(0, variance.y);
    p.z = center.z + vtkMath::Gaussian(0, variance.z);
    return p;
}


std::vector<Point3D> parseTrajectoryFile(const std::string& filePath, int n) {
    std::vector<Point3D> trajectory;
    std::ifstream file(filePath);
    std::string line;
    std::regex pattern(R"(x : -?\d+(\.\d+)?, y : -?\d+(\.\d+)?, z : -?\d+(\.\d+)?)");
    std::regex pattern2(R"([xyz:,\n])");
    std::smatch matches;

    int count = 0;
    while (std::getline(file, line) && (n <= 0 || count < n)) {
        if (std::regex_search(line, matches, pattern)) {
            Point3D point;
            std::istringstream iss(matches.str());
            char ignore;
            iss    >> ignore >> ignore >> point.x >>\
            ignore >> ignore >> ignore >> point.y >>\
            ignore >> ignore >> ignore >> point.z;
            trajectory.push_back(point);
            ++count;
        }
    }

    return trajectory;
}

void addVarianceToChart(vtkChartXYZ* chart, const std::vector<Point3D>& est_trajectory, const std::vector<Point3D>& variances, double p = 1.0, double s = 1.0) {
    vtkNew<vtkTable> table;
    vtkNew<vtkFloatArray> arrX, arrY, arrZ;
    vtkNew<vtkUnsignedCharArray> colors;
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    arrX->SetName("X");
    arrY->SetName("Y");
    arrZ->SetName("Z");

    double maxVariance = 0.0;
    for (const auto& var : variances) {
        double varianceMagnitude = sqrt(var.x * var.x + var.y * var.y + var.z * var.z);
        maxVariance = std::max(maxVariance, varianceMagnitude);
    }

    for (size_t i = 0; i < est_trajectory.size(); ++i) {
        if (vtkMath::Random() > s) continue;

        int numPoints = static_cast<int>(vtkMath::Random() * p);
        for (int j = 0; j < numPoints; ++j) {
            Point3D variancePoint = generateRandomPointAround(est_trajectory[i], variances[i]);

            arrX->InsertNextValue(variancePoint.x);
            arrY->InsertNextValue(variancePoint.y);
            arrZ->InsertNextValue(variancePoint.z);

            double normalizedVariance = sqrt(variances[i].x * variances[i].x + variances[i].y * variances[i].y + variances[i].z * variances[i].z) / maxVariance;
            double r = normalizedVariance * 255;
            double b = 255 * (1.0 - normalizedVariance);
            unsigned char color[3] =    {static_cast<unsigned char>(r),
                                        0,
                                        static_cast<unsigned char>(b)};
            colors->InsertNextTypedTuple(color);
        }
    }

    table->AddColumn(arrX);
    table->AddColumn(arrY);
    table->AddColumn(arrZ);
    table->GetRowData()->AddArray(colors);

    vtkNew<vtkPlotPoints3D> plot;
    plot->GetPen()->SetWidth(1.5);
    plot->SetInputData(table);

    chart->AddPlot(plot);
}



void addToChart(vtkChartXYZ* chart, const std::vector<Point3D>& points, const std::array<unsigned char, 3>& color) {
    vtkNew<vtkTable> table;
    vtkNew<vtkFloatArray> arrX, arrY, arrZ;
    arrX->SetName("X");
    arrY->SetName("Y");
    arrZ->SetName("Z");

    for (const auto& point : points) {
        arrX->InsertNextValue(point.x);
        arrY->InsertNextValue(point.y);
        arrZ->InsertNextValue(point.z);
    }

    table->AddColumn(arrX);
    table->AddColumn(arrY);
    table->AddColumn(arrZ);

    vtkNew<vtkPlotLine3D> plot;
    plot->SetInputData(table);

    plot->GetPen()->SetColorF(static_cast<double>(color[0]) / 255.0,
                              static_cast<double>(color[1]) / 255.0,
                              static_cast<double>(color[2]) / 255.0);

    plot->GetPen()->SetWidth(1.0);

    chart->AddPlot(plot);
}


class MouseWheelCallback : public vtkCommand {
public:
    static MouseWheelCallback *New() {
        return new MouseWheelCallback;
    }

    double increment = 0.05;

    void Execute(vtkObject *caller, unsigned long eventId, void *callData) override {
        if (eventId == vtkCommand::MouseWheelForwardEvent) {
            p += increment;
            s += increment;
        } else if (eventId == vtkCommand::MouseWheelBackwardEvent) {
            p = std::max(0.0, p - increment);
            s = std::max(0.0, s - increment);
        }
        replotVariance();
    }

    double p = increment;
    double s = increment;
    vtkChartXYZ* chart;
    std::vector<Point3D> truth;
    std::vector<Point3D> estimate;
    std::vector<Point3D> variance;

    void replotVariance() {
        chart->ClearPlots();
        addVarianceToChart(chart, estimate, variance, p, s);
        addToChart(chart, estimate, { 0, 255, 0 });
        addToChart(chart, truth, { 255, 0, 0 });

    }
};
 
int main(int, char*[])
{
    vtkObject::GlobalWarningDisplayOff();
    vtkNew<vtkChartXYZ> chart;
    vtkNew<vtkContextView> view;
    view->GetRenderWindow()->SetSize(600, 500);
    view->GetScene()->AddItem(chart);

    chart->SetMargins({ 40, 40, 40, 40 });

    chart->SetFitToScene(true);
    auto estimate = parseTrajectoryFile("predict.txt", -1);
    auto truth = parseTrajectoryFile("truth.txt", -1);
    auto variance = parseTrajectoryFile("variances.txt", -1);

    addToChart(chart, truth, { 255, 0, 0 });
    addToChart(chart, estimate, { 0, 255, 0 });
    addVarianceToChart(chart, estimate, variance);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    vtkNew<MouseWheelCallback> callback;
    callback->chart = chart.GetPointer();
    callback->estimate = estimate;
    callback->variance = variance;
    callback->truth = truth;
    view->GetInteractor()->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
    view->GetInteractor()->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);

    view->GetRenderWindow()->SetMultiSamples(0);
    view->GetInteractor()->Initialize();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    vtkContextMouseEvent mouseEvent;
    mouseEvent.SetInteractor(view->GetInteractor());
    vtkVector2f pos;
    vtkVector2f lastPos;

    // rotate
    mouseEvent.SetButton(vtkContextMouseEvent::LEFT_BUTTON);
    lastPos.Set(114, 55);
    mouseEvent.SetLastScenePos(lastPos);
    pos.Set(174, 121);
    mouseEvent.SetScenePos(pos);

    vtkVector2d sP(pos.Cast<double>().GetData());
    vtkVector2d lSP(lastPos.Cast<double>().GetData());

    vtkVector2d scenePos(mouseEvent.GetScenePos().Cast<double>().GetData());
    vtkVector2d lastScenePos(mouseEvent.GetLastScenePos().Cast<double>().GetData());
    chart->MouseMoveEvent(mouseEvent);

    // spin
    mouseEvent.SetButton(vtkContextMouseEvent::LEFT_BUTTON);
    mouseEvent.GetInteractor()->SetShiftKey(1);
    lastPos.Set(0, 0);
    mouseEvent.SetLastScenePos(lastPos);
    pos.Set(20, 10);
    mouseEvent.SetScenePos(pos);
    chart->MouseMoveEvent(mouseEvent);

    // zoom
    mouseEvent.SetButton(vtkContextMouseEvent::RIGHT_BUTTON);
    mouseEvent.GetInteractor()->SetShiftKey(0);
    lastPos.Set(0, 0);
    mouseEvent.SetLastScenePos(lastPos);
    pos.Set(0, 10);
    mouseEvent.SetScenePos(pos);
    chart->MouseMoveEvent(mouseEvent);

    // mouse wheel zoom
    chart->MouseWheelEvent(mouseEvent, -1);

    // pan
    mouseEvent.SetButton(vtkContextMouseEvent::LEFT_BUTTON);
    mouseEvent.GetInteractor()->SetShiftKey(1);
    lastPos.Set(0, 0);
    mouseEvent.SetLastScenePos(lastPos);
    pos.Set(100, 100);
    mouseEvent.SetScenePos(pos);
    chart->MouseMoveEvent(mouseEvent);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    view->GetRenderWindow()->Render();

    view->GetInteractor()->Start();

    return EXIT_SUCCESS;
}
