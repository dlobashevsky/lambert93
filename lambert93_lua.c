#include <stdlib.h>
#include <string.h>
#include <stdint.h>
//#include <time.h>

#include <lua.h>
#include <lauxlib.h>

#include "lambert93.h"

static int lua_latlon2xy(lua_State *L)
{
  if(lua_gettop(L)!=2 || !lua_isnumber(L,1) || !lua_isnumber(L,2))
  {
    fprintf(stderr,"%s failed, expected 2 numbers\n",__func__);
    lua_pushboolean(L,0);
    return 1;
  }
  double lat=lua_tonumber(L,1);
  double lon=lua_tonumber(L,2);
  double x=0;
  double y=0;
  lambert93_latlon2xy(lat, lon, &x, &y);
  lua_pushnumber(L,x);
  lua_pushnumber(L,y);
  return 2;
}

static int lua_xy2latlon(lua_State *L)
{
  if(lua_gettop(L)!=2 || !lua_isnumber(L,1) || !lua_isnumber(L,2))
  {
    fprintf(stderr,"%s failed, expected 2 numbers\n",__func__);
    lua_pushboolean(L,0);
    return 1;
  }
  double x=lua_tonumber(L,1);
  double y=lua_tonumber(L,2);
  double lat=0;
  double lon=0;
  int r=lambert93_xy2latlon(x, y, &lat, &lon);
  if(r)
  {
    fprintf(stderr,"%s failed, values are out of range\n",__func__);
    lua_pushboolean(L,0);
    return 1;
  }
  lua_pushnumber(L,lat);
  lua_pushnumber(L,lon);
  return 2;
}


static const struct luaL_Reg functions [] = {
    {"latlon2xy", lua_latlon2xy},
    {"xy2latlon", lua_xy2latlon},

    {NULL, NULL}
};

int luaopen_lambert93(lua_State *L)
{
    luaL_newlib(L,functions);
    return 1;
}

