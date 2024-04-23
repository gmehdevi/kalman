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
#include <array>
//regex
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>


struct Point3D {
    double x, y, z;
};

std::vector<Point3D> parseTrajectoryFile(const std::string& filePath, int n) {
    std::vector<Point3D> trajectory;
    std::ifstream file(filePath);
    std::string line;
    std::regex pattern(R"(x : -?\d+(\.\d+)?, y : -?\d+(\.\d+)?, z : -?\d+(\.\d+)?)");
    std::regex pattern2(R"([xyz:,\n])");
    std::smatch matches;

    int count = 0;
    while (std::getline(file, line) && (n <= 0 || count < n)) {
        std::cout << line << std::endl;
        if (std::regex_search(line, matches, pattern)) {
            Point3D point;
            std::istringstream iss(matches.str());
            std::cout << matches.str() << std::endl;
            char ignore;
            iss    >> ignore >> ignore >> point.x >>\
            ignore >> ignore >> ignore >> point.y >>\
            ignore >> ignore >> ignore >> point.z;
            trajectory.push_back(point);
            std::cout << point.x << " " << point.y << " " << point.z << std::endl;
            ++count;
        }
    }

    return trajectory;
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
    
    vtkNew<vtkUnsignedCharArray> colors;
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
    unsigned char vtkColor[3] = {color[0], color[1], color[2]};
    for (size_t i = 0; i < points.size(); ++i) {
        colors->InsertNextTypedTuple(vtkColor);
    }
    table->AddColumn(colors);
    plot->SetColors(colors);

    chart->AddPlot(plot);
}



int main(int, char*[])
{
    vtkNew<vtkChartXYZ> chart;
    vtkNew<vtkContextView> view;
    view->GetRenderWindow()->SetSize(600, 500);
    view->GetScene()->AddItem(chart);

    chart->SetMargins({ 40, 40, 40, 40 });

    chart->SetFitToScene(true);

    //   auto estimate = parseTrajectoryFile("predict.txt", 100);
    auto truth = parseTrajectoryFile("truth.txt", -1);
    //graph the points in the trajectory as a series of lines
    vtkNew<vtkTable> table3;
    for (auto& name : std::array<std::string, 3>{ "X", "Y", "Z" })
    {
        vtkNew<vtkFloatArray> arr;
        arr->SetName(name.c_str());
        table3->AddColumn(arr);
    }

    table3->SetNumberOfRows(truth.size());

    int pointindex = 0;
    for (auto& point : truth) {
        table3->SetValue(pointindex, 0, point.x);
        table3->SetValue(pointindex, 1, point.y);
        table3->SetValue(pointindex, 2, point.z);
        pointindex++;
    }

    vtkNew<vtkPlotLine3D> plot;
    plot->SetInputData(table3);
    chart->AddPlot(plot);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    view->GetRenderWindow()->SetMultiSamples(0);
    view->GetInteractor()->Initialize();
    view->GetRenderWindow()->Render();


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

    view->GetRenderWindow()->Render();

    view->GetInteractor()->Start();

    return EXIT_SUCCESS;
}
