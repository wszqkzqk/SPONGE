#pragma once

#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <highfive/highfive.hpp>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SpongeH5MD
{
struct ProtocolCVDefinition
{
    std::string name;
    std::string type;
    std::size_t dimension = 1;
    std::vector<std::pair<std::string, std::string>> runtime_parameters;
    std::vector<float> period;
    std::vector<float> sigma;
    std::vector<float> reference_coordinates;
};

class ProtocolCVH5Reader
{
   public:
    bool Open_Protocol(const std::string& file_path)
    {
        last_error_.clear();
        try
        {
            protocol_.reset(
                new HighFive::File(file_path, HighFive::File::ReadOnly));
            return true;
        }
        catch (const std::exception& error)
        {
            return Fail(std::string("failed to open protocol H5 file: ") +
                        error.what());
        }
    }

    bool Open_Restart(const std::string& file_path)
    {
        try
        {
            restart_.reset(
                new HighFive::File(file_path, HighFive::File::ReadOnly));
            return true;
        }
        catch (const std::exception& error)
        {
            return Fail(std::string("failed to open restart H5 file: ") +
                        error.what());
        }
    }

    bool Read_Definitions(std::size_t atom_count,
                          std::vector<ProtocolCVDefinition>* definitions)
    {
        if (definitions == nullptr)
        {
            return Fail("CV definition output pointer is null");
        }
        definitions->clear();
        if (protocol_ == nullptr) return Fail("protocol H5 reader is not open");
        if (!protocol_->exist("/cv")) return true;
        try
        {
            std::vector<std::string> names =
                protocol_->getGroup("/cv").listObjectNames();
            std::sort(names.begin(), names.end());
            for (const auto& name : names)
            {
                const std::string root = "/cv/" + name;
                if (name == "config" || protocol_->getObjectType(root) !=
                                            HighFive::ObjectType::Group)
                {
                    continue;
                }
                if (!Validate_Name(name, root) || !Read_Enabled(root)) continue;

                ProtocolCVDefinition definition;
                definition.name = name;
                definition.type = Read_Required_String(root, "type");
                definition.dimension = Read_Dimension(root);
                if (definition.dimension != 1)
                {
                    throw std::runtime_error(
                        root +
                        "/dimension must be 1 for the current scalar CV "
                        "runtime");
                }
                Add_Runtime_Parameter(&definition, "CV_type", definition.type,
                                      root + "/type");
                const auto atom_indices = Read_Atom_Indices(root, atom_count);
                if (!atom_indices.empty())
                {
                    Add_Runtime_Parameter(&definition, "atom",
                                          Join_Values(atom_indices),
                                          root + "/atom_indices");
                }
                Read_Parameter_Group(root, &definition);
                Read_Direct_Runtime_Fields(root, &definition);
                definition.period = Read_Optional_Float_Vector(
                    root + "/period", definition.dimension);
                definition.sigma = Read_Optional_Float_Vector(
                    root + "/sigma", definition.dimension);
                Read_Restart_Reference(root, atom_indices, &definition);
                Validate_Current_Runtime_Shape(atom_indices, definition);
                definitions->push_back(std::move(definition));
            }
            return true;
        }
        catch (const std::exception& error)
        {
            return Fail(std::string("failed to read typed CV definitions: ") +
                        error.what());
        }
    }

    const std::string& Last_Error() const { return last_error_; }

   private:
    bool Validate_Name(const std::string& name, const std::string& path)
    {
        if (name.empty() ||
            name.find_first_of(" \t\r\n,{}=/") != std::string::npos)
        {
            throw std::runtime_error(path + " has an invalid CV object name");
        }
        return true;
    }

    bool Read_Enabled(const std::string& root)
    {
        int enabled = 1;
        Read_Optional_Scalar(root, "enabled_default", &enabled);
        return enabled != 0;
    }

    std::string Read_Required_String(const std::string& root,
                                     const std::string& field)
    {
        std::string value;
        if (!Read_Optional_String(root, field, &value))
        {
            throw std::runtime_error("field is missing: " + root + "/" + field);
        }
        if (value.empty() || value.find_first_of("\r\n{}") != std::string::npos)
        {
            throw std::runtime_error(root + "/" + field +
                                     " contains an invalid value");
        }
        return value;
    }

    std::size_t Read_Dimension(const std::string& root)
    {
        long long dimension = 1;
        Read_Optional_Scalar(root, "dimension", &dimension);
        if (dimension <= 0 || static_cast<unsigned long long>(dimension) >
                                  std::numeric_limits<std::size_t>::max())
        {
            throw std::runtime_error(root +
                                     "/dimension must be a positive scalar");
        }
        return static_cast<std::size_t>(dimension);
    }

    std::vector<int> Read_Atom_Indices(const std::string& root,
                                       std::size_t atom_count)
    {
        const std::string path = root + "/atom_indices";
        if (!protocol_->exist(path))
        {
            std::string selection;
            if (Read_Optional_String(root, "selection_expression", &selection))
            {
                throw std::runtime_error(
                    root +
                    "/selection_expression is unresolved; SPONGE requires "
                    "atom_indices in sponge.input.v2");
            }
            return {};
        }
        const auto dims =
            protocol_->getDataSet(path).getSpace().getDimensions();
        if (dims.size() != 1 || dims[0] == 0)
        {
            throw std::runtime_error(path + " must have shape [n], n > 0");
        }
        std::vector<long long> raw(dims[0]);
        auto dataset = protocol_->getDataSet(path);
        if (H5Dread(dataset.getId(), H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, raw.data()) < 0)
        {
            throw std::runtime_error("failed to read " + path);
        }
        std::vector<int> result;
        result.reserve(raw.size());
        for (long long value : raw)
        {
            if (value < 0 ||
                value >
                    static_cast<long long>(std::numeric_limits<int>::max()) ||
                (atom_count != 0 &&
                 static_cast<std::size_t>(value) >= atom_count))
            {
                throw std::runtime_error(
                    path + " contains an out-of-range atom index");
            }
            result.push_back(static_cast<int>(value));
        }
        return result;
    }

    void Read_Parameter_Group(const std::string& root,
                              ProtocolCVDefinition* definition)
    {
        const std::string parameter_root = root + "/parameter";
        if (!protocol_->exist(parameter_root)) return;
        std::vector<std::string> fields =
            protocol_->getGroup(parameter_root).listObjectNames();
        std::sort(fields.begin(), fields.end());
        for (const auto& field : fields)
        {
            const std::string path = parameter_root + "/" + field;
            if (protocol_->getObjectType(path) != HighFive::ObjectType::Dataset)
            {
                throw std::runtime_error(path + " must be a typed dataset");
            }
            Validate_Parameter_Name(field, path);
            if (field == "CV_type" || field == "atom" ||
                field == "coordinate" || field == "period" ||
                field == "sigma" || field == "rotate" || field == "function" ||
                field == "min_padding" || field == "max_padding")
            {
                throw std::runtime_error(
                    path + " uses a reserved native CV field; write it at " +
                    root + "/" + field);
            }
            Add_Runtime_Parameter(definition, field, Read_Dataset_Tokens(path),
                                  path);
        }
    }

    void Read_Direct_Runtime_Fields(const std::string& root,
                                    ProtocolCVDefinition* definition)
    {
        std::string value;
        if (Read_Optional_String(root, "function", &value))
        {
            Add_Runtime_Parameter(definition, "function", value,
                                  root + "/function");
        }
        int rotate = 0;
        if (Read_Optional_Scalar(root, "rotate", &rotate))
        {
            Add_Runtime_Parameter(definition, "rotate", rotate != 0 ? "1" : "0",
                                  root + "/rotate");
        }
        double numeric = 0.0;
        for (const char* field : {"min_padding", "max_padding"})
        {
            if (Read_Optional_Scalar(root, field, &numeric))
            {
                Add_Runtime_Parameter(definition, field, Format_Float(numeric),
                                      root + "/" + field);
            }
        }
    }

    void Read_Restart_Reference(const std::string& root,
                                const std::vector<int>& atom_indices,
                                ProtocolCVDefinition* definition)
    {
        const std::string path = "/parameters/restart/references/cv/" +
                                 definition->name + "/coordinate";
        if (restart_ == nullptr || !restart_->exist(path))
        {
            if (definition->type == "rmsd" &&
                !Has_Runtime_Parameter(*definition, "coordinate"))
            {
                throw std::runtime_error(path +
                                         " is required for a native rmsd CV");
            }
            return;
        }
        if (definition->type != "rmsd")
        {
            throw std::runtime_error(path +
                                     " is only supported for rmsd CV objects");
        }
        const auto dims = restart_->getDataSet(path).getSpace().getDimensions();
        if (dims != std::vector<std::size_t>{atom_indices.size(), 3})
        {
            throw std::runtime_error(path +
                                     " must have shape [atom_indices,3]");
        }
        std::vector<float> values(atom_indices.size() * 3);
        auto dataset = restart_->getDataSet(path);
        if (H5Dread(dataset.getId(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, values.data()) < 0)
        {
            throw std::runtime_error("failed to read " + path);
        }
        Validate_Finite(values, path);
        definition->reference_coordinates = values;
        Add_Runtime_Parameter(definition, "coordinate", Join_Values(values),
                              path);
    }

    void Validate_Current_Runtime_Shape(const std::vector<int>& atom_indices,
                                        const ProtocolCVDefinition& definition)
    {
        std::size_t expected = 0;
        if (definition.type == "position_x" ||
            definition.type == "position_y" ||
            definition.type == "position_z" ||
            definition.type == "scaled_position_x" ||
            definition.type == "scaled_position_y" ||
            definition.type == "scaled_position_z")
            expected = 1;
        else if (definition.type == "distance" ||
                 definition.type == "displacement_x" ||
                 definition.type == "displacement_y" ||
                 definition.type == "displacement_z")
            expected = 2;
        else if (definition.type == "angle")
            expected = 3;
        else if (definition.type == "dihedral")
            expected = 4;
        else if (definition.type == "rmsd" && atom_indices.empty())
            throw std::runtime_error("native rmsd CV requires atom_indices");
        if (expected != 0 && atom_indices.size() != expected)
        {
            throw std::runtime_error("/cv/" + definition.name +
                                     "/atom_indices does not match CV type " +
                                     definition.type);
        }
    }

    std::vector<float> Read_Optional_Float_Vector(const std::string& path,
                                                  std::size_t count)
    {
        if (!protocol_->exist(path)) return {};
        const auto dims =
            protocol_->getDataSet(path).getSpace().getDimensions();
        if (dims != std::vector<std::size_t>{count})
        {
            throw std::runtime_error(path + " must match CV dimension");
        }
        std::vector<float> values(count);
        auto dataset = protocol_->getDataSet(path);
        if (H5Dread(dataset.getId(), H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                    H5P_DEFAULT, values.data()) < 0)
        {
            throw std::runtime_error("failed to read " + path);
        }
        Validate_Finite(values, path);
        return values;
    }

    std::string Read_Dataset_Tokens(const std::string& path)
    {
        auto dataset = protocol_->getDataSet(path);
        const auto dims = dataset.getSpace().getDimensions();
        std::size_t count = 1;
        for (std::size_t dim : dims) count *= dim;
        if (count == 0) throw std::runtime_error(path + " must not be empty");
        const H5T_class_t type_class =
            H5Tget_class(dataset.getDataType().getId());
        if (type_class == H5T_INTEGER)
        {
            std::vector<long long> values(count);
            if (H5Dread(dataset.getId(), H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, values.data()) < 0)
                throw std::runtime_error("failed to read " + path);
            return Join_Values(values);
        }
        if (type_class == H5T_FLOAT)
        {
            std::vector<double> values(count);
            if (H5Dread(dataset.getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, values.data()) < 0)
                throw std::runtime_error("failed to read " + path);
            Validate_Finite(values, path);
            return Join_Values(values);
        }
        if (type_class == H5T_STRING)
        {
            std::vector<std::string> values;
            if (dims.empty())
            {
                std::string value;
                dataset.read(value);
                values.push_back(value);
            }
            else if (dims.size() == 1)
            {
                dataset.read(values);
            }
            else
            {
                throw std::runtime_error(path +
                                         " string data must be scalar or 1D");
            }
            for (const auto& value : values)
            {
                if (value.empty() ||
                    value.find_first_of("\r\n{}") != std::string::npos)
                    throw std::runtime_error(path +
                                             " contains an invalid string");
            }
            return Join_Values(values);
        }
        throw std::runtime_error(path +
                                 " must use integer, float, or string data");
    }

    void Add_Runtime_Parameter(ProtocolCVDefinition* definition,
                               const std::string& name,
                               const std::string& value,
                               const std::string& source_path)
    {
        Validate_Parameter_Name(name, source_path);
        if (value.empty())
            throw std::runtime_error(source_path + " must not be empty");
        auto existing =
            std::find_if(definition->runtime_parameters.begin(),
                         definition->runtime_parameters.end(),
                         [&](const auto& item) { return item.first == name; });
        if (existing == definition->runtime_parameters.end())
        {
            definition->runtime_parameters.push_back({name, value});
        }
        else if (existing->second != value)
        {
            throw std::runtime_error(source_path +
                                     " conflicts with another CV field");
        }
    }

    bool Has_Runtime_Parameter(const ProtocolCVDefinition& definition,
                               const std::string& name) const
    {
        return std::any_of(definition.runtime_parameters.begin(),
                           definition.runtime_parameters.end(),
                           [&](const auto& item)
                           { return item.first == name; });
    }

    void Validate_Parameter_Name(const std::string& name,
                                 const std::string& path)
    {
        if (name.empty() ||
            name.find_first_of(" \t\r\n,{}=") != std::string::npos)
        {
            throw std::runtime_error(path +
                                     " has an invalid runtime parameter name");
        }
    }

    template <typename T>
    bool Read_Optional_Scalar(const std::string& root, const std::string& field,
                              T* value)
    {
        const std::string path = root + "/" + field;
        if (protocol_->exist(path))
        {
            protocol_->getDataSet(path).read(*value);
            return true;
        }
        auto group = protocol_->getGroup(root);
        if (group.hasAttribute(field))
        {
            group.getAttribute(field).read(*value);
            return true;
        }
        return false;
    }

    bool Read_Optional_String(const std::string& root, const std::string& field,
                              std::string* value)
    {
        return Read_Optional_Scalar(root, field, value);
    }

    template <typename T>
    static std::string Join_Values(const std::vector<T>& values)
    {
        std::ostringstream output;
        output << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            if (i != 0) output << ' ';
            output << values[i];
        }
        return output.str();
    }

    static std::string Format_Float(double value)
    {
        if (!std::isfinite(value))
            throw std::runtime_error("CV numeric field is not finite");
        std::ostringstream output;
        output << std::setprecision(std::numeric_limits<double>::max_digits10)
               << value;
        return output.str();
    }

    template <typename T>
    static void Validate_Finite(const std::vector<T>& values,
                                const std::string& path)
    {
        for (T value : values)
            if (!std::isfinite(value))
                throw std::runtime_error(path + " contains a non-finite value");
    }

    bool Fail(const std::string& message)
    {
        last_error_ = message;
        return false;
    }

    std::unique_ptr<HighFive::File> protocol_;
    std::unique_ptr<HighFive::File> restart_;
    std::string last_error_;
};
}  // namespace SpongeH5MD
