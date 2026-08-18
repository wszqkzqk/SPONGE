#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "third_party/toml/toml.h"

namespace
{
void Expect(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(1);
    }
}

void ExpectValue(const std::map<std::string, std::string>& commands,
                 const char* key, const char* expected)
{
    const auto iter = commands.find(key);
    Expect(iter != commands.end(), key);
    if (iter->second != expected)
    {
        std::cerr << "Expected " << key << " = " << expected << ", got "
                  << iter->second << '\n';
        std::exit(1);
    }
}
}  // namespace

int main()
{
    std::map<std::string, std::string> commands;
    std::string error;
    const std::string aliases = R"toml(
[write.interval]
information = 500
trajectory = 5000
mdout = 250
restart = 500000
)toml";

    Expect(sponge::toml_wrap::ParseAndFlatten(aliases, "aliases.in.spg.toml",
                                              &commands, &error),
           error.c_str());
    ExpectValue(commands, "write_information_interval", "500");
    ExpectValue(commands, "write_trajectory_interval", "5000");
    ExpectValue(commands, "write_mdout_interval", "250");
    ExpectValue(commands, "write_restart_file_interval", "500000");
    Expect(commands.count("write_interval_trajectory") == 0,
           "Alias key should not remain after normalization");

    const std::string conflict = R"toml(
[write]
trajectory_interval = 5000

[write.interval]
trajectory = 6000
)toml";
    commands.clear();
    error.clear();
    Expect(!sponge::toml_wrap::ParseAndFlatten(conflict, "conflict.in.spg.toml",
                                               &commands, &error),
           "Expected conflicting canonical and alias keys to fail");
    Expect(error.find("write_trajectory_interval") != std::string::npos,
           "Conflict error should name the canonical write key");

    const std::string mc_barostat = R"toml(
[barostat.monte_carlo]
update_interval = 100
initial_ratio = 0.002
accept_rate_low = 25
)toml";
    commands.clear();
    error.clear();
    Expect(sponge::toml_wrap::ParseAndFlatten(
               mc_barostat, "mc_barostat.in.spg.toml", &commands, &error),
           error.c_str());
    ExpectValue(commands, "monte_carlo_barostat_update_interval", "100");
    ExpectValue(commands, "monte_carlo_barostat_initial_ratio", "0.002");
    ExpectValue(commands, "monte_carlo_barostat_accept_rate_low", "25");
    Expect(
        commands.count("barostat_monte_carlo_update_interval") == 0,
        "Barostat Monte Carlo alias key should not remain after normalization");

    const std::string mc_conflict = R"toml(
[monte_carlo_barostat]
update_interval = 100

[barostat.monte_carlo]
update_interval = 200
)toml";
    commands.clear();
    error.clear();
    Expect(!sponge::toml_wrap::ParseAndFlatten(
               mc_conflict, "mc_conflict.in.spg.toml", &commands, &error),
           "Expected conflicting Monte Carlo barostat alias keys to fail");
    Expect(
        error.find("monte_carlo_barostat_update_interval") != std::string::npos,
        "Conflict error should name the canonical Monte Carlo barostat key");

    return 0;
}
