#ifndef EXCEPTIONHANDLER
#define EXCEPTIONHANDLER

#include <exception>
#include "gsl/gsl_errno.h"
#include "CppLib.h"

enum class ExceptionType
{
    FreezeIn,
    KineticDec_gammaLessThanH,
    KineticDec_gammaGreatThanH,
    GSL_error,
    Non_defined,
    WidthTooLarge
};
// type KineticDec_gammaLessThanH means that gamma < H and that the kinetic decoupling temperature is lower than the minimal value tested in the dichotomie for its estimation
// type KineticDec_gammaGreatThanH means that gamma > H and that the kinetic decoupling temperature is higher than the maximal value tested in the dichotomie for its estimation

class ExceptionHandler : public std::exception
{
public:
    ExceptionHandler() : _msg((char *)"WARNING << exception : "), _func(NULL), _line(0), _excType(ExceptionType::Non_defined) {};
    ExceptionHandler(char *msg) : _msg(msg), _func(NULL), _line(0), _excType(ExceptionType::Non_defined) {};
    ExceptionHandler(char *msg, int line, char *func) : _msg(msg), _func(func), _line(line), _excType(ExceptionType::Non_defined){};
    ExceptionHandler(char *msg, int line, char *func, ExceptionType exc) : _msg(msg), _func(func), _line(line), _excType(exc){};
    ExceptionHandler(char *msg, ExceptionType exc) : _msg(msg), _func(NULL), _line(0), _excType(exc){};
    ExceptionHandler(char *msg, int line, ExceptionType exc) : _msg(msg),  _func(NULL), _line(line), _excType(exc){};

    virtual const char *what() const throw() { return _msg; }
    const char *get_func() const { return _func; };
    int get_line() const { return _line; };
    ExceptionType get_exception_type() const { return _excType; };

private:
    char *_msg, *_func;
    int _line;
    ExceptionType _excType;
};

static void my_gsl_handler (const char * reason, const char * file,  int line, int gsl_errno)
{
    std::string reason_str = std::string(reason), file_str = std::string(file);
    std::stringstream line_str;
    line_str << line;
    std::string error_message =  "ERROR GSL : " + reason_str + " in " + file_str + " at line : " + line_str.str();
    
    if(NUM_THREADS == 1)
        std::cout << error_message <<  std::endl;

    char msg_array[error_message.size()];
    strcpy(msg_array, error_message.c_str());
    ExceptionHandler exc(msg_array, line, ExceptionType::GSL_error);

    throw exc;
}   


#endif
