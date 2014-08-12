#include <memory>
#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/particle.hpp"
#include "../examples/boris/particle_2d.hpp"

using namespace ::testing;

typedef ParticleEncapsulation<double, double> ParticleT;

typedef Position2DEncapsulation<double, double> Position2D;
typedef Velocity2DEncapsulation<double, double> Velocity2D;
typedef Acceleration2DEncapsulation<double, double> Acceleration2D;
typedef Particle2DEncapsulation<double, double> Particle2DT;


TEST(BorisParticle, Instantiation)
{
  ParticleT default_ctor;
  EXPECT_THAT(default_ctor.DIM, 0);
  EXPECT_THAT(default_ctor.charge, DoubleEq(1.0));
  EXPECT_THAT(default_ctor.mass, DoubleEq(1.0));

  ParticleT special_ctor = ParticleT(0.5, 1.0);
  EXPECT_THAT(special_ctor.mass, DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge, DoubleEq(1.0));
}

/*
 * Position2DEncapsulation
 */
TEST(BorisParticlePosition2D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Position2D>::value);
  EXPECT_TRUE(is_destructible<Position2D>::value);
  EXPECT_TRUE(is_copy_constructible<Position2D>::value);
  EXPECT_TRUE(is_move_constructible<Position2D>::value);

  Position2D default_pos;
  EXPECT_THAT(default_pos.DIM, 2);
  EXPECT_THAT(default_pos.x, DoubleEq(0.0));
  EXPECT_THAT(default_pos.y, DoubleEq(0.0));

  shared_ptr<const Position2D> pos2 = make_shared<const Position2D>(0.5, 1.0);
  EXPECT_THAT(pos2->x, DoubleEq(0.5));
  EXPECT_THAT(pos2->y, DoubleEq(1.0));
}

TEST(BorisParticlePosition2D, Copyable)
{
  EXPECT_TRUE(is_copy_assignable<Position2D>::value);
  EXPECT_TRUE(is_move_assignable<Position2D>::value);

  shared_ptr<const Position2D> pos_original = make_shared<const Position2D>(0.5, 1.0);

  shared_ptr<Position2D> pos_copy1 = make_shared<Position2D>();
  pos_copy1->copy(pos_original);
  EXPECT_THAT(pos_copy1->x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy1->y, DoubleEq(1.0));

  Position2D pos_copy2;
  pos_copy2 = pos_original;
  EXPECT_THAT(pos_copy2.x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy2.y, DoubleEq(1.0));

  Position2D pos_copy3;
  pos_copy3 = *(pos_original.get());
  EXPECT_THAT(pos_copy3.x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy3.y, DoubleEq(1.0));
}

TEST(BorisParticlePosition2D, axpy)
{
  shared_ptr<Position2D> pos_2 = make_shared<Position2D>(0.5, 1.0);
  shared_ptr<const Position2D> other_pos = make_shared<const Position2D>(0.5, 1.0);

  pos_2->saxpy(1.0, other_pos);
  EXPECT_THAT(pos_2->x, DoubleEq(1.0));
  EXPECT_THAT(pos_2->y, DoubleEq(2.0));

  shared_ptr<Position2D> pos_1 = make_shared<Position2D>(0.5, 1.0);
  pos_1->saxpy(0.0, other_pos);
  EXPECT_THAT(pos_1->x, DoubleEq(0.5));
  EXPECT_THAT(pos_1->y, DoubleEq(1.0));

  shared_ptr<Position2D> pos_0 = make_shared<Position2D>(0.5, 1.0);
  pos_0->saxpy(-1.0, other_pos);
  EXPECT_THAT(pos_0->x, DoubleEq(0.0));
  EXPECT_THAT(pos_0->y, DoubleEq(0.0));
}

TEST(BorisParticlePosition2D, operators)
{
  Position2D pos1(0.1, 0.2);
  Position2D pos2(-0.1, -0.2);

  long double one = 1.0;

  Position2D twice = pos1 * 2.0;
  EXPECT_THAT(twice.x, DoubleEq(0.2));
  EXPECT_THAT(twice.y, DoubleEq(0.4));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position2D twice2 = 2.0 * pos1;
  EXPECT_THAT(twice2.x, DoubleEq(0.2));
  EXPECT_THAT(twice2.y, DoubleEq(0.4));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position2D add1 = pos1 + pos2;
  EXPECT_THAT(add1.x, DoubleEq(0.0));
  EXPECT_THAT(add1.y, DoubleEq(0.0));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));
  EXPECT_THAT(pos2.x, DoubleEq(-0.1));
  EXPECT_THAT(pos2.y, DoubleEq(-0.2));

  Position2D add2 = pos1 + one;
  EXPECT_THAT(add2.x, DoubleEq(1.1));
  EXPECT_THAT(add2.y, DoubleEq(1.2));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position2D add3 = one + pos1;
  EXPECT_THAT(add3.x, DoubleEq(1.1));
  EXPECT_THAT(add3.y, DoubleEq(1.2));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));
}


/*
 * Velocity2DEncapsulation
 */
TEST(BorisParticleVelocity2D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Velocity2D>::value);
  EXPECT_TRUE(is_copy_constructible<Velocity2D>::value);
  EXPECT_TRUE(is_move_constructible<Velocity2D>::value);
  EXPECT_TRUE(is_destructible<Velocity2D>::value);

  Velocity2D default_vel;
  EXPECT_THAT(default_vel.DIM, 2);
  EXPECT_THAT(default_vel.u, DoubleEq(0.0));
  EXPECT_THAT(default_vel.v, DoubleEq(0.0));

  shared_ptr<const Velocity2D> vel2 = make_shared<const Velocity2D>(0.5, 1.0);
  EXPECT_THAT(vel2->u, DoubleEq(0.5));
  EXPECT_THAT(vel2->v, DoubleEq(1.0));
}

TEST(BorisParticleVelocity2D, Copyable)
{
  EXPECT_TRUE(is_copy_assignable<Velocity2D>::value);
  EXPECT_TRUE(is_move_assignable<Velocity2D>::value);

  shared_ptr<Velocity2D> vel_copy = make_shared<Velocity2D>();
  shared_ptr<const Velocity2D> vel2 = make_shared<const Velocity2D>(0.5, 1.0);
  vel_copy->copy(vel2);
  EXPECT_THAT(vel_copy->u, DoubleEq(0.5));
  EXPECT_THAT(vel_copy->v, DoubleEq(1.0));
}

TEST(BorisParticleVelocity2D, axpy)
{
  shared_ptr<Velocity2D> vel_2 = make_shared<Velocity2D>(0.5, 1.0);
  shared_ptr<const Velocity2D> vel2 = make_shared<const Velocity2D>(0.5, 1.0);
  vel_2->saxpy(1.0, vel2);
  EXPECT_THAT(vel_2->u, DoubleEq(1.0));
  EXPECT_THAT(vel_2->v, DoubleEq(2.0));

  shared_ptr<Velocity2D> vel_1 = make_shared<Velocity2D>(0.5, 1.0);
  vel_1->saxpy(0.0, vel2);
  EXPECT_THAT(vel_1->u, DoubleEq(0.5));
  EXPECT_THAT(vel_1->v, DoubleEq(1.0));

  shared_ptr<Velocity2D> vel_0 = make_shared<Velocity2D>(0.5, 1.0);
  vel_0->saxpy(-1.0, vel2);
  EXPECT_THAT(vel_0->u, DoubleEq(0.0));
  EXPECT_THAT(vel_0->v, DoubleEq(0.0));
}

TEST(BorisParticleVelocity2D, convert)
{
  Velocity2D vel1(0.1, 0.2);

  Position2D converted = vel1.convert(dt<double>(2.0));
  EXPECT_THAT(converted.x, DoubleEq(0.2));
  EXPECT_THAT(converted.y, DoubleEq(0.4));
}

TEST(BorisParticleVelocity2D, operators)
{
  Velocity2D vel1(0.1, 0.2);
  Velocity2D vel2(-0.1, -0.2);

  long double one = 1.0;

  Velocity2D twice = vel1 * 2.0;
  EXPECT_THAT(twice.u, DoubleEq(0.2));
  EXPECT_THAT(twice.v, DoubleEq(0.4));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity2D twice2 = 2.0 * vel1;
  EXPECT_THAT(twice2.u, DoubleEq(0.2));
  EXPECT_THAT(twice2.v, DoubleEq(0.4));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity2D add1 = vel1 + vel2;
  EXPECT_THAT(add1.u, DoubleEq(0.0));
  EXPECT_THAT(add1.v, DoubleEq(0.0));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));
  EXPECT_THAT(vel2.u, DoubleEq(-0.1));
  EXPECT_THAT(vel2.v, DoubleEq(-0.2));

  Velocity2D add2 = vel1 + one;
  EXPECT_THAT(add2.u, DoubleEq(1.1));
  EXPECT_THAT(add2.v, DoubleEq(1.2));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity2D add3 = one + vel1;
  EXPECT_THAT(add3.u, DoubleEq(1.1));
  EXPECT_THAT(add3.v, DoubleEq(1.2));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));
}


/*
 * Acceleration2DEncapsulation
 */
TEST(BorisParticleAcceleration2D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Acceleration2D>::value);
  EXPECT_TRUE(is_copy_constructible<Acceleration2D>::value);
  EXPECT_TRUE(is_move_constructible<Acceleration2D>::value);
  EXPECT_TRUE(is_destructible<Acceleration2D>::value);

  Acceleration2D default_accel;
  EXPECT_THAT(default_accel.DIM, 2);
  EXPECT_THAT(default_accel.a, DoubleEq(0.0));
  EXPECT_THAT(default_accel.b, DoubleEq(0.0));

  shared_ptr<const Acceleration2D> accel2 = make_shared<const Acceleration2D>(0.5, 1.0);
  EXPECT_THAT(accel2->a, DoubleEq(0.5));
  EXPECT_THAT(accel2->b, DoubleEq(1.0));
}

TEST(BorisParticleAcceleration2D, Copyable)
{
  shared_ptr<Acceleration2D> accel_copy = make_shared<Acceleration2D>();
  shared_ptr<const Acceleration2D> accel2 = make_shared<const Acceleration2D>(0.5, 1.0);
  accel_copy->copy(accel2);
  EXPECT_THAT(accel_copy->a, DoubleEq(0.5));
  EXPECT_THAT(accel_copy->b, DoubleEq(1.0));
}

TEST(BorisParticleAcceleration2D, axpy)
{
  shared_ptr<Acceleration2D> accel_2 = make_shared<Acceleration2D>(0.5, 1.0);
  shared_ptr<const Acceleration2D> accel2 = make_shared<const Acceleration2D>(0.5, 1.0);
  accel_2->saxpy(1.0, accel2);
  EXPECT_THAT(accel_2->a, DoubleEq(1.0));
  EXPECT_THAT(accel_2->b, DoubleEq(2.0));

  shared_ptr<Acceleration2D> accel_1 = make_shared<Acceleration2D>(0.5, 1.0);
  accel_1->saxpy(0.0, accel2);
  EXPECT_THAT(accel_1->a, DoubleEq(0.5));
  EXPECT_THAT(accel_1->b, DoubleEq(1.0));

  shared_ptr<Acceleration2D> accel_0 = make_shared<Acceleration2D>(0.5, 1.0);
  accel_0->saxpy(-1.0, accel2);
  EXPECT_THAT(accel_0->a, DoubleEq(0.0));
  EXPECT_THAT(accel_0->b, DoubleEq(0.0));
}

TEST(BorisParticleAcceleration2D, convert)
{
  Acceleration2D accel1(0.1, 0.2);

  Position2D converted_pos = accel1.convert(dtdt<double>(2.0));
  EXPECT_THAT(converted_pos.x, DoubleEq(0.2));
  EXPECT_THAT(converted_pos.y, DoubleEq(0.4));

  Velocity2D converted_vel = accel1.convert(dt<double>(2.0));
  EXPECT_THAT(converted_vel.u, DoubleEq(0.2));
  EXPECT_THAT(converted_vel.v, DoubleEq(0.4));
}

TEST(BorisParticleAcceleration2D, operators)
{
  Acceleration2D accel1(0.1, 0.2);
  Acceleration2D accel2(-0.1, -0.2);

  long double one = 1.0;

  Acceleration2D twice = accel1 * 2.0;
  EXPECT_THAT(twice.a, DoubleEq(0.2));
  EXPECT_THAT(twice.b, DoubleEq(0.4));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration2D twice2 = 2.0 * accel1;
  EXPECT_THAT(twice2.a, DoubleEq(0.2));
  EXPECT_THAT(twice2.b, DoubleEq(0.4));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration2D add1 = accel1 + accel2;
  EXPECT_THAT(add1.a, DoubleEq(0.0));
  EXPECT_THAT(add1.b, DoubleEq(0.0));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));
  EXPECT_THAT(accel2.a, DoubleEq(-0.1));
  EXPECT_THAT(accel2.b, DoubleEq(-0.2));

  Acceleration2D add2 = accel1 + one;
  EXPECT_THAT(add2.a, DoubleEq(1.1));
  EXPECT_THAT(add2.b, DoubleEq(1.2));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration2D add3 = one + accel1;
  EXPECT_THAT(add3.a, DoubleEq(1.1));
  EXPECT_THAT(add3.b, DoubleEq(1.2));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));
}


/*
 * Particle2DEncapsulation
 */
TEST(BorisParticle2D, Instantiation)
{
  EXPECT_TRUE(is_default_constructible<Particle2DT>::value);
  EXPECT_TRUE(is_copy_constructible<Particle2DT>::value);
  EXPECT_TRUE(is_move_constructible<Particle2DT>::value);
  EXPECT_TRUE(is_destructible<Particle2DT>::value);

  Particle2DT default_ctor;
  EXPECT_THAT(default_ctor.DIM, 2);
  EXPECT_THAT(default_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.mass, DoubleEq(1.0));
  EXPECT_THAT(default_ctor.charge, DoubleEq(1.0));

  Particle2DT special_ctor = Particle2DT(0.5, 1.0);
  EXPECT_THAT(special_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.mass, DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge, DoubleEq(1.0));
}

TEST(BorisParticle2D, Copyable)
{
  shared_ptr<Particle2DT> original_ptr = make_shared<Particle2DT>(0.1, 0.2);
  original_ptr->pos->x = 0.3;
  EXPECT_THAT(original_ptr->pos->x, DoubleEq(0.3));
  original_ptr->pos->y = 0.4;
  EXPECT_THAT(original_ptr->pos->y, DoubleEq(0.4));

  original_ptr->vel->u = 0.5;
  EXPECT_THAT(original_ptr->vel->u, DoubleEq(0.5));
  original_ptr->vel->v = 0.6;
  EXPECT_THAT(original_ptr->vel->v, DoubleEq(0.6));

  original_ptr->accel->a = 0.7;
  EXPECT_THAT(original_ptr->accel->a, DoubleEq(0.7));
  original_ptr->accel->b = 0.8;
  EXPECT_THAT(original_ptr->accel->b, DoubleEq(0.8));
  Particle2DT expl_copy;

  expl_copy.copy(const_pointer_cast<const Particle2DT>(original_ptr));
  EXPECT_THAT(expl_copy.pos->x, DoubleEq(0.3));
  EXPECT_THAT(expl_copy.pos->y, DoubleEq(0.4));
  EXPECT_THAT(expl_copy.vel->u, DoubleEq(0.5));
  EXPECT_THAT(expl_copy.vel->v, DoubleEq(0.6));
  EXPECT_THAT(expl_copy.accel->a, DoubleEq(0.7));
  EXPECT_THAT(expl_copy.accel->b, DoubleEq(0.8));
  EXPECT_THAT(expl_copy.mass, DoubleEq(0.1));
  EXPECT_THAT(expl_copy.charge, DoubleEq(0.2));

  Particle2DT original = Particle2DT(0.1, 0.2);
  original.pos->x = 0.3;
  original.pos->y = 0.4;
  original.vel->u = 0.5;
  original.vel->v = 0.6;
  original.accel->a = 0.7;
  original.accel->b = 0.8;
  Particle2DT impl_copy = original;
  EXPECT_THAT(impl_copy.pos->x, DoubleEq(0.3));
  EXPECT_THAT(impl_copy.pos->y, DoubleEq(0.4));
  EXPECT_THAT(impl_copy.vel->u, DoubleEq(0.5));
  EXPECT_THAT(impl_copy.vel->v, DoubleEq(0.6));
  EXPECT_THAT(impl_copy.accel->a, DoubleEq(0.7));
  EXPECT_THAT(impl_copy.accel->b, DoubleEq(0.8));
  EXPECT_THAT(impl_copy.mass, DoubleEq(0.1));
  EXPECT_THAT(impl_copy.charge, DoubleEq(0.2));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
