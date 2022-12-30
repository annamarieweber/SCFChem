#ifndef TIMED_FUNCTIONAL
#define TIMED_FUNCTIONAL

#include <iostream>
#include <chrono>
#include <functional>
#include <armadillo>
#include <typeinfo>
#include "cndo2.h"
#include "cluster.h"

namespace util
{
    namespace timing
    {
        /**
         * @brief A class that wraps a function or function pointer and measures the elapsed time when the function is called.
         *
         * @tparam F the functional to time.
         * @tparam Args The argument types of the function or function pointer.
         */
        template <typename F, typename... Args>
        class TimedFunctional

        {
        public:
            F fn; ///< The function or function pointer to be wrapped.
            string fn_name;

            TimedFunctional(F &&val);

            TimedFunctional(F &&val, string name);

            /**
             * @brief Calls the wrapped function or function pointer and measures the elapsed time.
             *
             * @param args The arguments to be passed to the wrapped function or function pointer.
             * @return The return value of the wrapped function or function pointer.
             */
            typename std::invoke_result<F, Args...>::type operator()(Args... args);

            /**
             * @brief Prints the elapsed time and the name of the wrapped function or function pointer to the console.
             *
             * @param elapsed The elapsed time to be printed.
             */
            void printTiming(std::chrono::duration<double> elapsed);
        };

        template <typename F, typename... Args>
        TimedFunctional<F, Args...>::TimedFunctional(F &&val)
        {
            fn = val;
            fn_name = typeid(fn).name();
        }

        template <typename F, typename... Args>
        TimedFunctional<F, Args...>::TimedFunctional(F &&val, string name)
        {
            fn = val;
            fn_name = name;
        }

        template <typename F, typename... Args>
        typename std::invoke_result<F, Args...>::type TimedFunctional<F, Args...>::operator()(Args... args)
        {
            std::chrono::time_point<std::chrono::system_clock> start, end;
            std::chrono::duration<double> elapsed;

            // starting clock
            start = std::chrono::system_clock::now();
            // calling function and saving result
            auto result = fn(args...);
            // stopping clock
            end = std::chrono::system_clock::now();
            // calculating elapsed time
            elapsed = end - start;
            // print timing output
            TimedFunctional::printTiming(elapsed);
            return result;
        }

        template <typename F, typename... Args>
        void TimedFunctional<F, Args...>::printTiming(std::chrono::duration<double> elapsed)
        {
            std::cout << "Executed Call  " << fn_name << ": " << elapsed.count() << " secs" << std::endl;
        }
    }
}

#endif