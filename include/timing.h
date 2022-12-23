#ifndef TIMER_UTIL
#define TIMER_UTIL

#include <iostream>
#include <chrono>
#include <functional>
#include <typeinfo>

namespace util
{
    namespace timing
    {
        /**
         * @brief A struct that wraps a function or function pointer and measures the elapsed time when the function is called.
         *
         * @tparam returnType The return type of the function or function pointer.
         * @tparam Args The argument types of the function or function pointer.
         */
        template <typename returnType, typename... Args>
        struct TimedFunctional
        {
            std::function<returnType(Args...)> fn; ///< The function or function pointer to be wrapped.

            /**
             * @brief Calls the wrapped function or function pointer and measures the elapsed time.
             *
             * @param args The arguments to be passed to the wrapped function or function pointer.
             * @return The return value of the wrapped function or function pointer.
             */
            returnType operator()(Args... args);

            /**
             * @brief Prints the elapsed time and the name of the wrapped function or function pointer to the console.
             *
             * @param elapsed The elapsed time to be printed.
             */
            void printTiming(std::chrono::duration<double> elapsed);
        };
    }
}
#endif // TIMER_UTIL
