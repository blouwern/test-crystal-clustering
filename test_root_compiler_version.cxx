#include <iostream>
#include <print>

auto test_root_compiler_version() -> int{
    std::println("Your root current compiler version: {}", __cplusplus);
    // std::printf("Your root current compiler version: %ld\n", __cplusplus);
    return 0;
}
