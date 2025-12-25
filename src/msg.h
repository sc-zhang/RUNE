//
// Created by zsc on 2023/11/6.
//

#ifndef RUNE_MSG_H
#define RUNE_MSG_H
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

class msg {
private:
  void print(const std::string &str, const std::string &color_string) const {
    time_t now = time(nullptr);
    tm *tm_t = localtime(&now);
    if (this->use_cout) {
      std::cout << color_string << std::setw(2) << std::setfill('0')
                << tm_t->tm_hour << ":" << std::setw(2) << std::setfill('0')
                << tm_t->tm_min << ":" << std::setw(2) << std::setfill('0')
                << tm_t->tm_sec << "] \033[0m " << str << "\n";
    } else {
      std::cerr << color_string << std::setw(2) << std::setfill('0')
                << tm_t->tm_hour << ":" << std::setw(2) << std::setfill('0')
                << tm_t->tm_min << ":" << std::setw(2) << std::setfill('0')
                << tm_t->tm_sec << "] \033[0m " << str << "\n";
    }
  };
  bool use_cout;

public:
  explicit msg(bool use_cout = true) { this->use_cout = use_cout; }
  void info(const std::string &str) { this->print(str, "\033[32m ["); };
  void warn(const std::string &str) { this->print(str, "\033[33m ["); };
  void err(const std::string &str) { this->print(str, "\033[31m ["); };
};

#endif // RUNE_MSG_H
