layout(local_size_x_id = 0, local_size_y_id = 1, local_size_z_id = 2) in;

struct Vec3
{
    float x, y, z;
};

struct LTMat
{
    float a11, a21, a22, a31, a32, a33;
};

struct Boundary
{
    uint policy;
    float a11, a21, a22, a31, a32, a33;
    float r11, r21, r22, r31, r32, r33;
};

Vec3 v3(float x, float y, float z) { return Vec3(x, y, z); }
Vec3 v3_add(Vec3 a, Vec3 b) { return Vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
Vec3 v3_sub(Vec3 a, Vec3 b) { return Vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
Vec3 v3_scale(Vec3 a, float s) { return Vec3(a.x * s, a.y * s, a.z * s); }
Vec3 v3_mul(Vec3 a, Vec3 b) { return Vec3(a.x * b.x, a.y * b.y, a.z * b.z); }
Vec3 v3_neg(Vec3 a) { return Vec3(-a.x, -a.y, -a.z); }
float v3_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
float v3_len(Vec3 a) { return sqrt(v3_dot(a, a)); }
Vec3 v3_floor(Vec3 a) { return Vec3(floor(a.x), floor(a.y), floor(a.z)); }
Vec3 v3_cross(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x);
}

float copysignf(float mag, float sgn)
{
    return uintBitsToFloat((floatBitsToUint(mag) & 0x7FFFFFFFu) |
                           (floatBitsToUint(sgn) & 0x80000000u));
}

Vec3 v3_mul_lt(Vec3 v, float a11, float a21, float a22, float a31, float a32,
               float a33)
{
    return Vec3(v.x * a11 + v.y * a21 + v.z * a31, v.y * a22 + v.z * a32,
                v.z * a33);
}

Vec3 get_displacement(Vec3 a, Vec3 b, Boundary bd)
{
    Vec3 dr = v3_sub(a, b);
    if (bd.policy == 0u) return dr;
    Vec3 s = v3_mul_lt(dr, bd.r11, bd.r21, bd.r22, bd.r31, bd.r32, bd.r33);
    s = v3_floor(v3_add(s, Vec3(0.5, 0.5, 0.5)));
    return v3_sub(dr, v3_mul_lt(s, bd.a11, bd.a21, bd.a22, bd.a31, bd.a32,
                                bd.a33));
}

Vec3 get_displacement_periodic(Vec3 a, Vec3 b, Boundary bd)
{
    Vec3 dr = v3_sub(a, b);
    Vec3 s = v3_mul_lt(dr, bd.r11, bd.r21, bd.r22, bd.r31, bd.r32, bd.r33);
    s = v3_floor(v3_add(s, Vec3(0.5, 0.5, 0.5)));
    return v3_sub(dr, v3_mul_lt(s, bd.a11, bd.a21, bd.a22, bd.a31, bd.a32,
                                bd.a33));
}

struct AtomGroup
{
    int atom_numbers;
    int ghost_numbers;
    uint serial_lo;
    uint serial_hi;
};

Vec3 wrap_coordinate(Vec3 c, Boundary bd)
{
    if (bd.policy == 0u) return c;
    Vec3 s = v3_mul_lt(c, bd.r11, bd.r21, bd.r22, bd.r31, bd.r32, bd.r33);
    s = v3_floor(s);
    return v3_sub(c, v3_mul_lt(s, bd.a11, bd.a21, bd.a22, bd.a31, bd.a32,
                               bd.a33));
}

LTMat virial_from_force_dis(Vec3 f, Vec3 dr)
{
    LTMat m;
    m.a11 = f.x * dr.x;
    m.a21 = f.x * dr.y + f.y * dr.x;
    m.a22 = f.y * dr.y;
    m.a31 = f.x * dr.z + f.z * dr.x;
    m.a32 = f.y * dr.z + f.z * dr.y;
    m.a33 = f.z * dr.z;
    return m;
}

#define afadd(slot, v) atomicAdd(slot, v)
#define aadd_v3(slot, v)                       \
    do                                         \
    {                                          \
        atomicAdd((slot).x, (v).x);            \
        atomicAdd((slot).y, (v).y);            \
        atomicAdd((slot).z, (v).z);            \
    } while (false)
#define aadd_lt(slot, m)            \
    do                              \
    {                               \
        atomicAdd((slot).a11, (m).a11); \
        atomicAdd((slot).a21, (m).a21); \
        atomicAdd((slot).a22, (m).a22); \
        atomicAdd((slot).a31, (m).a31); \
        atomicAdd((slot).a32, (m).a32); \
        atomicAdd((slot).a33, (m).a33); \
    } while (false)

uvec2 u64_add(uvec2 a, uvec2 b)
{
    uint lo = a.x + b.x;
    return uvec2(lo, a.y + b.y + (lo < a.x ? 1u : 0u));
}
uvec2 u64_shr2(uvec2 a) { return uvec2((a.x >> 2) | (a.y << 30), a.y >> 2); }

uvec4 philox_round(uvec4 counter, uvec2 key)
{
    uint lo0, hi0, lo1, hi1;
    umulExtended(0xD2511F53u, counter.x, hi0, lo0);
    umulExtended(0xCD9E8D57u, counter.z, hi1, lo1);
    return uvec4(hi1 ^ counter.y ^ key.x, lo1, hi0 ^ counter.w ^ key.y, lo0);
}

uvec4 philox_block(uvec2 block_offset, uvec2 stream, uvec2 seed)
{
    uvec4 counter = uvec4(block_offset.x, block_offset.y, stream.x, stream.y);
    uvec2 key = seed;
    for (int round = 0; round < 10; round++)
    {
        counter = philox_round(counter, key);
        if (round != 9) key += uvec2(0x9E3779B9u, 0xBB67AE85u);
    }
    return counter;
}

vec4 philox_normal4(uvec2 scalar_offset, uvec2 stream, uvec2 seed)
{
    uvec4 values = philox_block(u64_shr2(scalar_offset), stream, seed);
    uint substate = scalar_offset.x & 3u;
    uint words[4];
    for (uint index = 0u; index < 4u; index++)
    {
        uint absolute = substate + index;
        if (absolute < 4u)
        {
            words[index] = values[absolute];
        }
        else
        {
            uvec4 next = philox_block(
                u64_add(u64_shr2(scalar_offset), uvec2(1u, 0u)), stream, seed);
            words[index] = next[absolute - 4u];
        }
    }
    const float scale = 1.0 / 4294967295.0;
    float u1 = float(words[0]) * scale;
    float u2 = float(words[1]) * scale;
    float u3 = float(words[2]) * scale;
    float u4 = float(words[3]) * scale;
    u1 = u1 > 0.0 ? sqrt(-2.0 * log(u1)) : 0.0;
    u2 = 6.283185308 * u2;
    u3 = u3 > 0.0 ? sqrt(-2.0 * log(u3)) : 0.0;
    u4 = 6.283185308 * u4;
    return vec4(u1 * cos(u2), u1 * sin(u2), u3 * cos(u4), u3 * sin(u4));
}

uint philox_uint32(uvec2 scalar_offset, uvec2 stream, uvec2 seed)
{
    uvec4 values = philox_block(u64_shr2(scalar_offset), stream, seed);
    return values[scalar_offset.x & 3u];
}

shared float sg_reduce[32];
shared int sg_reduce_i[32];

float block_sum(float v)
{
    uint lane = gl_SubgroupInvocationID;
    uint sid = gl_SubgroupID;
    v = subgroupAdd(v);
    if (lane == 0u) sg_reduce[sid] = v;
    barrier();
    uint total = gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z;
    uint nsub = (total + gl_SubgroupSize - 1u) / gl_SubgroupSize;
    if (sid == 0u)
    {
        v = lane < nsub ? sg_reduce[lane] : 0.0;
        v = subgroupAdd(v);
        if (lane == 0u) sg_reduce[0] = v;
    }
    barrier();
    return sg_reduce[0];
}

int block_sum_int(int v)
{
    uint lane = gl_SubgroupInvocationID;
    uint sid = gl_SubgroupID;
    v = subgroupAdd(v);
    if (lane == 0u) sg_reduce_i[sid] = v;
    barrier();
    uint total = gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z;
    uint nsub = (total + gl_SubgroupSize - 1u) / gl_SubgroupSize;
    if (sid == 0u)
    {
        v = lane < nsub ? sg_reduce_i[lane] : 0;
        v = subgroupAdd(v);
        if (lane == 0u) sg_reduce_i[0] = v;
    }
    barrier();
    return sg_reduce_i[0];
}

struct Kahan
{
    float sum, c;
};
void kahan_add(inout Kahan k, float v)
{
    float y = v - k.c;
    float t = k.sum + y;
    k.c = (t - k.sum) - y;
    k.sum = t;
}

LTMat lt_scale(LTMat m, float s)
{
    return LTMat(m.a11 * s, m.a21 * s, m.a22 * s, m.a31 * s, m.a32 * s,
                 m.a33 * s);
}

LTMat lt_sub(LTMat a, LTMat b)
{
    return LTMat(a.a11 - b.a11, a.a21 - b.a21, a.a22 - b.a22, a.a31 - b.a31,
                 a.a32 - b.a32, a.a33 - b.a33);
}

int get_lj_type(int a, int b)
{
    int hi = max(a, b);
    int lo = min(a, b);
    return (hi * (hi + 1) >> 1) + lo;
}

// Abramowitz-Stegun 7.1.26, |epsilon| <= 1.2e-7
float erfcf(float x)
{
    float ax = abs(x);
    float t = 1.0 / (1.0 + 0.5 * ax);
    float tau =
        t * exp(-ax * ax - 1.26551223 +
                t * (1.00002368 +
                     t * (0.37409196 +
                          t * (0.09678418 +
                               t * (-0.18628806 +
                                    t * (0.27886807 +
                                         t * (-1.13520398 +
                                              t * (1.48851587 +
                                                   t * (-0.82215223 +
                                                        t * 0.17087277)))))))));
    return x >= 0.0 ? tau : 2.0 - tau;
}

#define warp_afadd(slot, v)                                      \
    do                                                           \
    {                                                            \
        float wsum_ = subgroupAdd(v);                            \
        if (gl_SubgroupInvocationID == 0u) atomicAdd(slot, wsum_); \
    } while (false)
#define warp_aadd_v3(slot, v)                                \
    do                                                       \
    {                                                        \
        Vec3 wsum_ = v;                                      \
        wsum_.x = subgroupAdd(wsum_.x);                      \
        wsum_.y = subgroupAdd(wsum_.y);                      \
        wsum_.z = subgroupAdd(wsum_.z);                      \
        if (gl_SubgroupInvocationID == 0u) aadd_v3(slot, wsum_); \
    } while (false)
#define warp_aadd_lt(slot, v)                                \
    do                                                       \
    {                                                        \
        LTMat wsum_ = v;                                     \
        wsum_.a11 = subgroupAdd(wsum_.a11);                  \
        wsum_.a21 = subgroupAdd(wsum_.a21);                  \
        wsum_.a22 = subgroupAdd(wsum_.a22);                  \
        wsum_.a31 = subgroupAdd(wsum_.a31);                  \
        wsum_.a32 = subgroupAdd(wsum_.a32);                  \
        wsum_.a33 = subgroupAdd(wsum_.a33);                  \
        if (gl_SubgroupInvocationID == 0u) aadd_lt(slot, wsum_); \
    } while (false)
