#include <Random123/threefry.h>
#include <Random123/u01fixedpt.h>

float rand(float16* neutron, uint tid) {
  float counter = (*neutron).sb;

  threefry4x32_key_t k = {{tid, 0xdecafbad, 0xfacebead, 0x12345678}};
  threefry4x32_ctr_t c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};

  union {
    threefry4x32_ctr_t c;
    int4 i;
  } u;

  c.v[0]+=counter;
  u.c = threefry4x32(c, k);

  (*neutron).sb+=1;

  return (float)u01fixedpt_closed_closed_32_24(u.i.x);
}