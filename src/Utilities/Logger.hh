#ifndef SPHERAL_LOGGER_HH
#define SPHERAL_LOGGER_HH

#include "config.hh"

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cstdio>

// Check if the debug macro is enabled by CMake
#ifdef SPHERAL_ENABLE_LOGGER


// DebugLogStream class to enable << operator for logging
class DebugLogStream {
public:
    // Constructor captures file, line, and function information
    DebugLogStream() {
        // Start the stream with the standard prefix
        m_stream << "[DEBUG] ";
    }

    // Destructor flushes the accumulated message to stderr
    ~DebugLogStream() {
        std::fprintf(stderr, "%s\n", m_stream.str().c_str());
    }

    // Overload the << operator to allow streaming various types
    template<typename T>
    DebugLogStream& operator<<(const T& msg) {
        m_stream << msg;
        return *this;
    }

    // Special overload for manipulators like std::endl (though not strictly needed
    // as we flush in destructor, it makes the syntax feel more natural if used)
    DebugLogStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        manip(m_stream); // Apply the manipulator to the internal stream
        return *this;
    }

private:
    std::stringstream m_stream;
};

// Macro to create a temporary DebugLogStream object
// This object will print its content when it goes out of scope (end of statement)
#define DEBUG_LOG DebugLogStream()

#else // SPHERAL_ENABLE_LOGGER not defined

// If debug is not enabled, the macro expands to a null stream,
// effectively removing all debug logging from the compiled code.
// This uses a dummy class to consume the stream operators without doing anything.
class NullDebugLogStream {
public:
    // Consume any type
    template<typename T>
    NullDebugLogStream& operator<<(const T& msg) {
        (void)msg; // Suppress unused variable warning
        return *this;
    }
    // Consume manipulators
    NullDebugLogStream& operator<<(std::ostream& (*manip)(std::ostream&)) {
        (void)manip; // Suppress unused variable warning
        return *this;
    }
};

// Define DEBUG_LOG to return a temporary NullDebugLogStream object
#define DEBUG_LOG NullDebugLogStream()

#endif // SPHERAL_ENABLE_LOGGER


#endif // SPHERAL_LOGGER_HH
