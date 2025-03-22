//
// Created by zsc on 2023/11/7.
//

#ifndef RUNE_RUNE_H
#define RUNE_RUNE_H
namespace rune {
class FLAG {
public:
  static const uint32_t UNKNOWN = 0;
};
class MASK {
public:
  static const uint64_t ID_MASK = 0xFFFFFFFF00000000;
  static const uint64_t POS_MASK = 0x00000000FFFFFFFF;
};
} // namespace rune
#endif // RUNE_RUNE_H
