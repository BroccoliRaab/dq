#ifndef DQ_H
#define DQ_H

#ifndef DQ_SINF
#include <math.h>
#define DQ_SINF sinf
#endif

#ifndef DQ_COSF
#include <math.h>
#define DQ_COSF cosf
#endif

#ifndef DQ_SQRTF
#include <math.h>
#define DQ_SQRTF sqrtf
#endif

typedef struct dq_quat_s dq_quat_t;
typedef struct dq_vec3_s dq_vec3_t;
typedef struct dq_mat3_s dq_mat3_t;
typedef struct dq_dual_quat_s dq_dual_quat_t;

struct dq_quat_s {
    float w;
    float x;
    float y;
    float z;
};

struct dq_dual_quat_s {
    float rw;
    float rx;
    float ry;
    float rz;
    float dw;
    float dx;
    float dy;
    float dz;
};

struct dq_vec3_s
{
    float x;
    float y;
    float z;
};

struct dq_mat3_s
{
    float ax, ay, az;
    float bx, by, bz;
    float cx, cy, cz;
};
static const dq_quat_t dq_quat_unit = {
    1.0, 0.0, 0.0, 0.0
};

static const dq_dual_quat_t dq_dual_quat_unit = {
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static inline float 
dq_vec3_mag(
    dq_vec3_t v
)
{
    return DQ_SQRTF(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline dq_vec3_t 
dq_vec3_normalize(
    const dq_vec3_t v
)
{
    const float imag = 1.0/dq_vec3_mag(v);
    return (dq_vec3_t){v.x*imag, v.y*imag, v.z*imag};
}

static inline float 
dq_vec3_dot(
    dq_vec3_t a,
    dq_vec3_t b
){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline dq_vec3_t 
dq_vec3_cross(
    dq_vec3_t a,
    dq_vec3_t b
){
    return (dq_vec3_t){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

static inline dq_quat_t 
dq_quat_scalar(
    float s
)
{
    return (dq_quat_t){s, 0.0, 0.0, 0.0};
}

static inline dq_quat_t 
dq_quat_negate(
    dq_quat_t q
)
{
    return (dq_quat_t){-q.w, -q.x, -q.y, -q.z};
}

static inline dq_quat_t 
dq_quat_conjugate(
    dq_quat_t q
)
{
    return (dq_quat_t){q.x, -q.y, -q.z, -q.w};
}


static inline float 
dq_quat_mag2(
    dq_quat_t q
)
{
    return q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
}

static inline float 
dq_quat_mag(
    dq_quat_t q
)
{
    return DQ_SQRTF(q.x) + DQ_SQRTF(q.y) + DQ_SQRTF(q.z) + DQ_SQRTF(q.w);
}

static inline dq_quat_t 
dq_quat_smul(
    dq_quat_t q,
    float s
)
{
    return (dq_quat_t){s*q.w, s*q.x, s*q.y, s*q.z};
}

static inline dq_quat_t 
dq_quat_divs(
    dq_quat_t q,
    float s
)
{
    return (dq_quat_t){q.w/s, q.x/s, q.y/s, q.z/s};
}

static inline dq_quat_t
dq_quat_normalize(
    dq_quat_t q
)
{
    return dq_quat_divs(
        q,
        dq_quat_mag(q)
    );
}

static inline dq_quat_t
dq_quat_inverse(
    dq_quat_t q
)
{
    return dq_quat_divs(
        dq_quat_conjugate(q),
        dq_quat_mag(q)
    );
}

static inline float
dq_quat_dot(
    dq_quat_t a,
    dq_quat_t b
)
{
    return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w;
}

static inline dq_quat_t
dq_quat_cross(
    dq_quat_t a,
    dq_quat_t b
){

    const float dot = dq_vec3_dot(
        (dq_vec3_t){a.y, a.z, a.w},
        (dq_vec3_t){b.y, b.z, b.w}
    );
    const dq_vec3_t cross = dq_vec3_cross(
        (dq_vec3_t){a.y, a.z, a.w},
        (dq_vec3_t){b.y, b.z, b.w}
    );

    return (dq_quat_t){
        a.x * b.x - dot,
        a.x * b.y + b.x * a.y + cross.x,
        a.x * b.z + b.x * a.z + cross.y,
        a.x * b.w + b.x * a.w + cross.z
    };
    

}

static inline dq_quat_t
dq_quat_mul(
    dq_quat_t a,
    dq_quat_t b
)
{
    return (dq_quat_t){
        (b.w * a.w) - (b.x * a.x) - (b.y * a.y) - (b.z * a.z),
        (b.w * a.x) + (b.x * a.w) + (b.y * a.z) - (b.z * a.y),
        (b.w * a.y) + (b.y * a.w) + (b.z * a.x) - (b.x * a.z),
        (b.w * a.z) + (b.z * a.w) + (b.x * a.y) - (b.y * a.x)
    };
}
static inline dq_quat_t
dq_quat_add(
    dq_quat_t a,
    dq_quat_t b
)
{
    return (dq_quat_t){
        (b.w + a.w),
        (b.x + a.x), 
        (b.y + a.y),
        (b.z + a.z) 
    };
}

static inline dq_quat_t 
dq_quat_from_axis_angle(
    float angle,
    dq_vec3_t axis
)
{
    const dq_vec3_t axis_norm = dq_vec3_normalize(axis);
    const float axis_coef = DQ_SINF(angle * 0.5);
    return (dq_quat_t){
        DQ_COSF(angle * 0.5),
        axis.x * axis_coef,
        axis.y * axis_coef,
        axis.z * axis_coef
    };
}

static inline dq_mat3_t
dq_quat_to_mat3(
    dq_quat_t q
)
{
    const float s = q.x;
    const float x = q.y;
    const float y = q.z;
    const float z = q.w;
    const float x2 = 2.0 * x * x;
    const float y2 = 2.0 * y * z;
    const float z2 = 2.0 * z * z;
    const float xy = 2.0 * x * y;
    const float yz = 2.0 * z * y;
    const float xz = 2.0 * z * x;
    const float sz = 2.0 * z * s;
    const float sx = 2.0 * x * s;
    const float sy = 2.0 * y * s;
    return (dq_mat3_t){
        1.0 - y2 - z2, xy - sz, xz + sy,
        xy + sz, 1.0 - x2 - z2, yz - sx, 
        xz - sy, yz - sx, 1.0 - x2 - y2
    };
}

static inline dq_dual_quat_t
dq_dual_quat_mul(
    dq_dual_quat_t a,
    dq_dual_quat_t b
)
{
    dq_dual_quat_t res;
    ((dq_quat_t *) &res)[0] = dq_quat_mul(
        ((dq_quat_t *) &a)[0],
        ((dq_quat_t *) &b)[0]
    );
    ((dq_quat_t *) &res)[0] = dq_quat_add(
        dq_quat_mul(
            ((dq_quat_t *) &a)[0],
            ((dq_quat_t *) &b)[1]
        ),
        dq_quat_mul(
            ((dq_quat_t *) &a)[1],
            ((dq_quat_t *) &b)[0]
        )
    );
    return res;
}

static inline dq_dual_quat_t
dq_dual_quat_pose(
    dq_quat_t rot,
    dq_vec3_t pos
)
{
    dq_dual_quat_t pose;
    ((dq_quat_t *) &pose)[0] = rot;
    ((dq_quat_t *) &pose)[1] = dq_quat_divs(
        dq_quat_mul(
            (dq_quat_t){0.0, pos.x, pos.y, pos.z},
            rot
        ), 
        2.0
    );
    return pose;
}

static inline dq_vec3_t
dq_dual_quat_tranform(
    dq_dual_quat_t dq,
    dq_vec3_t point
){
    dq_dual_quat_t coeff;
    ((dq_quat_t*) &coeff)[0] = dq_quat_conjugate(((dq_quat_t *) &dq)[0]);
    ((dq_quat_t*) &coeff)[1] = dq_quat_negate(
        dq_quat_conjugate(((dq_quat_t *) &dq)[0])
    );

    const dq_dual_quat_t res = dq_dual_quat_mul(
        dq_dual_quat_mul(
            dq,
            (dq_dual_quat_t){1.0, 0.0, 0.0, 0.0, 0.0, point.x, point.y, point.z}
        ),
        coeff
    );
    return (dq_vec3_t){res.ry, res.rz, res.rw};
}


#endif /* DQ_H */
