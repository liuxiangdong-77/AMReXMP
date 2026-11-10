#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <vtkDataArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>

namespace {

struct Options {
  std::filesystem::path directory;
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double dt{1.0};
  std::filesystem::path output{"particle_stress_history.csv"};
};

struct Frame {
  long long step{0};
  std::filesystem::path path;
};

std::optional<long long> frame_step(const std::filesystem::path& path) {
  const std::string stem = path.stem().string();
  const std::string prefix = "particle_data";
  if (stem.rfind(prefix, 0) != 0 || stem.size() == prefix.size())
    return std::nullopt;
  const std::string suffix = stem.substr(prefix.size());
  if (!std::all_of(suffix.begin(), suffix.end(),
                   [](unsigned char value) { return std::isdigit(value); }))
    return std::nullopt;
  return std::stoll(suffix);
}

bool parse_options(int argc, char** argv, Options* options) {
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    auto next_value = [&](const std::string& name) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << name << '\n';
        return nullptr;
      }
      return argv[++i];
    };

    if (argument == "--dir") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->directory = value;
    } else if (argument == "--x") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->x = std::stod(value);
    } else if (argument == "--y") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->y = std::stod(value);
    } else if (argument == "--z") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->z = std::stod(value);
    } else if (argument == "--dt") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->dt = std::stod(value);
    } else if (argument == "--output") {
      const char* value = next_value(argument);
      if (value == nullptr) return false;
      options->output = value;
    } else if (argument == "--help" || argument == "-h") {
      return false;
    } else {
      std::cerr << "Unknown argument: " << argument << '\n';
      return false;
    }
  }
  return !options->directory.empty() && options->dt > 0.0;
}

void print_usage(const char* executable) {
  std::cout
      << "Usage:\n  " << executable
      << " --dir VTP_DIRECTORY --x X --y Y --z Z [--dt TIME_STEP]"
         " [--output FILE.csv]\n\n"
      << "Example:\n  " << executable
      << " --dir results-amr-two-level-sine/amr-two-level-sine-000"
         " --x 5.125 --y 5.125 --z 5.125 --dt 1e-5"
         " --output excitation_stress.csv\n";
}

}  // namespace

int main(int argc, char** argv) {
  Options options;
  if (!parse_options(argc, argv, &options)) {
    print_usage(argv[0]);
    return 1;
  }

  std::vector<Frame> frames;
  if (!std::filesystem::is_directory(options.directory)) {
    std::cerr << "Directory does not exist: " << options.directory << '\n';
    return 1;
  }
  for (const auto& entry : std::filesystem::directory_iterator(options.directory)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".vtp")
      continue;
    const auto step = frame_step(entry.path());
    if (step.has_value()) frames.push_back(Frame{*step, entry.path()});
  }
  std::sort(frames.begin(), frames.end(),
            [](const Frame& lhs, const Frame& rhs) {
              return lhs.step < rhs.step;
            });
  if (frames.empty()) {
    std::cerr << "No particle_data*.vtp files were found\n";
    return 1;
  }

  std::ofstream output(options.output);
  if (!output) {
    std::cerr << "Cannot create output file: " << options.output << '\n';
    return 1;
  }
  output << std::setprecision(16)
         << "step,time,sxx,syy,szz,sxy,syz,szx,von_mises,x,y,z\n";

  for (const auto& frame : frames) {
    vtkNew<vtkXMLPolyDataReader> reader;
    reader->SetFileName(frame.path.string().c_str());
    reader->Update();

    auto* polydata = reader->GetOutput();
    auto* points = polydata->GetPoints();
    auto* stresses = polydata->GetPointData()->GetArray("stresses");
    if (points == nullptr || stresses == nullptr ||
        stresses->GetNumberOfComponents() < 9) {
      std::cerr << "Missing points or nine-component stresses in "
                << frame.path << '\n';
      return 1;
    }

    double nearest_distance = std::numeric_limits<double>::max();
    vtkIdType nearest_id = -1;
    double nearest_point[3]{};
    for (vtkIdType id = 0; id < points->GetNumberOfPoints(); ++id) {
      double point[3]{};
      points->GetPoint(id, point);
      const double dx = point[0] - options.x;
      const double dy = point[1] - options.y;
      const double dz = point[2] - options.z;
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (distance < nearest_distance) {
        nearest_distance = distance;
        nearest_id = id;
        std::copy(point, point + 3, nearest_point);
      }
    }

    double components[9]{};
    stresses->GetTuple(nearest_id, components);
    const double sxx = components[0];
    const double syy = components[4];
    const double szz = components[8];
    const double sxy = components[1];
    const double syz = components[5];
    const double szx = components[2];
    const double von_mises =
        std::sqrt(0.5 * ((sxx - syy) * (sxx - syy) +
                         (syy - szz) * (syy - szz) +
                         (szz - sxx) * (szz - sxx)) +
                  3.0 * (sxy * sxy + syz * syz + szx * szx));

    output << frame.step << ',' << frame.step * options.dt << ',' << sxx << ','
           << syy << ',' << szz << ',' << sxy << ',' << syz << ',' << szx
           << ',' << von_mises << ',' << nearest_point[0] << ','
           << nearest_point[1] << ',' << nearest_point[2] << '\n';
  }

  std::cout << "Wrote " << frames.size() << " frames to " << options.output
            << '\n';
  return 0;
}
