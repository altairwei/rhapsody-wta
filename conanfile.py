from conans import ConanFile, CMake, tools


class RhapsodyWtaConan(ConanFile):
    name = "rhapsody-wta"
    version = "0.0.1"
    license = "MIT"
    author = "Altair Wei altair_wei@outlook.com"
    url = "https://github.com/altairwei/rhapsody-wta.git"
    description = "Some tools for handling the output of Rhapsody WTA pipeline"
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_find_package"
    requires = "bamtools/2.5.1@altairwei/testing"

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin")
        self.copy("*.dylib*", dst="bin", src="lib")

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
