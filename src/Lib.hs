module Lib (
    Point(..),
    xgcd,
    modinv,
    multiply,
    divide,
    tangent,
    identify,
    dot,
    add,
    double,
    pointMultiply
) where

import Data.Bits
import Data.Bool


-- a simple data class to store (x, y) pairs
data Point = Point {
    x :: Integer,
    y :: Integer
}

-- the extended euclidean algorithm
xgcd :: Integer -> Integer -> (Integer, Integer, Integer)
xgcd 0 b = (b, 0, 1)
xgcd a b = (x1, x2, x3) where
    (gc, x, y) = xgcd (mod b a) a
    x1 = gc
    x2 = y - (b `div` a) * x
    x3 = x

-- inverse mod
modinv :: Integer -> Integer -> Integer
modinv a b
    | g /= 1 = error "g not equal to 0"
    | otherwise = (mod x b)
    where
        (g, x, y) = xgcd a b

-- multiply on a finite field
multiply :: Integer -> Integer -> Integer -> Integer
multiply a b p = mod (a * b) p

-- divide on a finite field
divide :: Integer -> Integer -> Integer -> Integer
divide a b p = (multiply m ec p)
    where
        ec = modinv b p
        m = mod a p

-- get the tangent of a point on a finite field
tangent :: (Point, Integer, Integer) -> Integer
tangent (point, a, p) = divide ((x point * x point * 3) + a) (y point * 2) p

-- get the identity of a point by prime
identify :: Integer -> Point
identify p = Point p 0

-- dots the curve
dot :: (Point, Point, Integer, Integer) -> Point
dot (p1, p2, m, p) = Point t q
    where
        v = mod (y p1 + p - mod (m * x p1) p) p
        t = mod (m * m + p - x p1 + p - x p2) p
        q = mod (p - (mod (m * t) p) + p - v) p

-- doubles a point
double :: (Point, Integer, Integer) -> Point
double (point, a, p)
    | x point == p = point
    | otherwise = dot (point, point, tan, p)
        where
            tan = tangent(point, a, p)

-- adds two points
add :: (Point, Point, Integer, Integer) -> Point
add (p1, p2, p, a)
    | x p1 == p = p2
    | x p2 == p = p1
    | (x p1 == x p2) && (y p1 == y p2) = double(p1, a, p)
    | (x p1 == x p2) = identify(p)
    | otherwise = dot (p1, p2, m, p)
        where
            m = divide (y p1 + p - y p2) (x p1 + p - x p2) p

-- r is identity of generator, q is generator, scalar multiplication
pointMultiply :: (Point, Integer, Integer, Point, Integer) -> Point
pointMultiply (r, 0, p, q, a) = r
pointMultiply (r, n, p, q, a)
    | ((.&.)  n 1 /= 0) = pointMultiply((add (r, q, p, a)), nm, p, double(q, a, p), a)
    | otherwise = pointMultiply(r, nm, p, double(q, a, p), a)
    where
        nm = shiftR n 1

-- signs data with private ket
sign :: (Integer, Integer, Point, Integer, Integer, Integer) -> Point
sign (pri, d, g, p, a, n) = Point r s
    where
        k = mod d n
        p1 = pointMultiply ((identify p), k, p, g, a)
        r = x p1 -- todo check not equal to zero
        s = (modinv k n) * mod (d + mod (pri * r) n) n -- todo check not equal to zero

verify :: (Point, Integer, Point, Point, Integer, Integer, Integer) -> Bool
verify (pub, d, sig, g, n, p, a)
    | r < 1 || r > n - 1 = False
    | s < 1 || s > n - 1 = False
    | otherwise = r == v
    where
        r = x sig
        s = y sig
        c = modinv s n
        u1 = mod (d * c) n
        u2 = mod (r * c) n
        xy = add ((pointMultiply ((identify p), u1, p, g, a)), (pointMultiply ((identify p), u1, p, pub, a)), p, a)
        v = mod (x xy) n

