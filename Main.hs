{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FunctionalDependencies #-}

import Data.IORef
data Dvar a = Dvar { value :: a, backprop :: (a, Int) -> IO a }

instance Num a => Num (Dvar a) where
    (+) (Dvar x dx) (Dvar y dy) = Dvar (x + y) (\(dz, j) -> do
        dx_result <- dx (dz, j)
        dy_result <- dy (dz, j)
        return (dx_result + dy_result))
    (-) (Dvar x dx) (Dvar y dy) = Dvar (x - y) (\(dz, j) -> do
        dx_result <- dx (dz, j)
        dy_result <- dy (-dz, j)
        return (dx_result + dy_result))
    (*) (Dvar x dx) (Dvar y dy) = Dvar (x * y) (\(dz, j) -> do
        dx_result <- dx (y * dz, j)
        dy_result <- dy (x * dz, j)
        return (dx_result + dy_result))
    negate (Dvar x dx) = Dvar (negate x) (\(dz, j) -> do
        result <- dx (negate dz, j)
        return result)
    abs (Dvar x dx) = Dvar (abs x) (\(dz, j) -> dx (signum x * dz, j))
    signum (Dvar x dx) = Dvar (signum x) (\(_, _) -> return 0)
    fromInteger n = Dvar (fromInteger n) (\(_, _) -> return 0)

instance Fractional a => Fractional (Dvar a) where
    (/) (Dvar x dx) (Dvar y dy) = Dvar (x / y) (\(dz, j) -> do
        dx_result <- dx (dz / y, j)
        dy_result <- dy ((-x / (y * y)) * dz, j)
        return (dx_result + dy_result))
    fromRational r = Dvar (fromRational r) (\(_, _) -> return 0)

instance Floating a => Floating (Dvar a) where
    (**) (Dvar x dx) (Dvar y dy) = Dvar (x ** y) (\(dz, j) -> do
        dx_result <- dx (y * (x ** (y - 1)) * dz, j)
        dy_result <- dy ((x ** y) * log x * dz, j)
        return (dx_result + dy_result))
    pi = Dvar pi (\(_, _) -> return 0)
    exp (Dvar x dx) = Dvar (exp x) (\(dz, j) -> do
        result <- dx (exp x * dz, j)
        return result)
    log (Dvar x dx) = Dvar (log x) (\(dz, j) -> do
        result <- dx (dz / x, j)
        return result)
    sin (Dvar x dx) = Dvar (sin x) (\(dz, j) -> do
        result <- dx (cos x * dz, j)
        return result)
    cos (Dvar x dx) = Dvar (cos x) (\(dz, j) -> do
        result <- dx (-(sin x) * dz, j)
        return result)
    asin (Dvar x dx) = Dvar (asin x) (\(dz, j) -> do
        result <- dx (dz / sqrt(1 - x*x), j)
        return result)
    acos (Dvar x dx) = Dvar (acos x) (\(dz, j) -> do
        result <- dx (-dz / sqrt(1 - x*x), j)
        return result)
    atan (Dvar x dx) = Dvar (atan x) (\(dz, j) -> do
        result <- dx (dz / (1 + x*x), j)
        return result)
    sinh (Dvar x dx) = Dvar (sinh x) (\(dz, j) -> do
        result <- dx (cosh x * dz, j)
        return result)
    cosh (Dvar x dx) = Dvar (cosh x) (\(dz, j) -> do
        result <- dx (sinh x * dz, j)
        return result)
    asinh (Dvar x dx) = Dvar (asinh x) (\(dz, j) -> do
        result <- dx (dz / sqrt(x*x + 1), j)
        return result)
    acosh (Dvar x dx) = Dvar (acosh x) (\(dz, j) -> do
        result <- dx (dz / sqrt(x*x - 1), j)
        return result)
    atanh (Dvar x dx) = Dvar (atanh x) (\(dz, j) -> do
        result <- dx (dz / (1 - x*x), j)
        return result)
    
jac :: (Num a, ToList b a) => ([Dvar a] -> b) -> [a] -> IO [[a]]
jac f inputs = do
    let n = length inputs 
        m = length (toList (f (map (\x -> Dvar x (\(_, _) -> return 0)) inputs)))
    jacobianRef <- newIORef (replicate m (replicate n 0))
    let dvars = zipWith (lift jacobianRef) [0..] inputs
        results = toList (f dvars)
    mapM_ (\(j, result) -> backprop result (1, j)) (zip [0..] results)
    readIORef jacobianRef

lift :: Num a => IORef [[a]] -> Int -> a -> Dvar a
lift jacobianRef i x = Dvar x (\(dz, j) -> do
    modifyIORef jacobianRef (\jac -> 
        take j jac ++ 
        [take i (jac !! j) ++ [dz + ((jac !! j) !! i)] ++ drop (i+1) (jac !! j)] ++ 
        drop (j+1) jac)
    return dz)

class ToList a b | a -> b where
    toList :: a -> [Dvar b]

instance ToList (Dvar a) a where
    toList x = [x]

instance ToList [Dvar a] a where
    toList = id

f :: Floating a => [a] -> a
f [x, y, z] = x*x + x*y + sin z

f2 [x] = 2*x**3
f3 [x,y] = x+y
f4 [x,y] = f3 [x,y] + f2 [x]

c x = 2*x**3
c1 x y = x+y
c2 x y = c1 x y + c x

g [x,y] = y*sin x
g1 [x,y] = (g [x,y])**2
g2 [x,y] = (g1 [x,y]) + y

h [x] = sin x
h1 [x,y] = h [x] + y
h2 [x,y] = (h1 [x,y]) * y

m x y = x*y
m1 x y = (m x y) + cos y
m2 [x,y] = x*(m1 x y)

n :: Floating a => [a] -> [a]
n [x, y] = [x*x + y, sin x + y*y]

main :: IO ()
main = do
    gradient <- jac f [2.0, 3.0, 1.0]
    print gradient
    gradient <- jac f2 [2.0]
    print gradient
    gradient <- jac f4 [5,1]
    print gradient
    gradient <- jac g2 [2,3]
    print gradient
    gradient <- jac h2 [2,3]
    print gradient
    gradient <- jac m2 [5,2]
    print gradient
    jacobian <- jac n [2.0, 3.0]
    print jacobian