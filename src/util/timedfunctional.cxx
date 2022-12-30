#include <iostream>
#include <chrono>
#include <functional>
#include <typeinfo>
#include "timedfunctional.h"

template <typename F, typename... Args>
class TimedFunctional
{
public:
    F fn; ///< The function or function pointer to be wrapped.

    /**
     * @brief Construct a new Timed Functional< F,  Args...>:: Timed Functional object
     *
     * @param val the function to wrap
     */
    TimedFunctional(F &&val)
    {
        fn = val;
        fn_name = typeid(fn).name();
    }

    /**
     * @brief Construct a new Timed Functional< F,  Args...>:: Timed Functional object
     *
     * @param val the function to wrap
     * @param name the name of the function
     */
    TimedFunctional(F &&val, string name);
    {
        fn = val;
        fn_name = name;
    }

    /**
     * @brief Calls the wrapped function or function pointer and measures the elapsed time.
     *
     * @param args The arguments to be passed to the wrapped function or function pointer.
     * @return The return value of the wrapped function or function pointer.
     */
    typename std::invoke_result<F, Args...>::type operator()(Args... args)
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

    /**
     * @brief Prints the elapsed time and the name of the wrapped function or function pointer to the console.
     *
     * @param elapsed The elapsed time to be printed.
     */
    void printTiming(std::chrono::duration<double> elapsed)
    {
        std::cout << "Executed Call  " << fn_name << ": " << elapsed.count() << " secs" << std::endl;
    }
};