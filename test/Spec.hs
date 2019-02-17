import Test.HUnit
import Test.Framework
import Test.Framework.Providers.HUnit
import Data.Monoid
import Control.Monad
import Lib


xgcdTest :: Assertion
xgcdTest = do
    assertEqual "XGCD case 1" (30, 1, -1) (xgcd 180 150)
    assertEqual "XGCD case 2" (1, -26, 31) (xgcd 180 151)
    assertEqual "XGCD case 3" (30, 1, -1) (xgcd 180 150)

modinvTest :: Assertion
modinvTest = do
    assertEqual "MODINV case 1" (9) (modinv 3 26)
    assertEqual "MODINV case 2" (19) (modinv 323 26)
    assertEqual "MODINV case 3" (57) (modinv 323 263)
    assertEqual "MODINV case 4" (18824) (modinv 323 26321)

multiplyTest :: Assertion
multiplyTest = do
    assertEqual "MUTPILY case 1" 6 (multiply 4 5 7)
    assertEqual "MUTPILY case 2" 20 (multiply 4 5 99)
    assertEqual "MUTPILY case 3" 60 (multiply 40 51 99)

divideTest :: Assertion
divideTest = do
    assertEqual "DIVIDE case 1" (17) (divide 40 1 23)
    assertEqual "DIVIDE case 2" (6) (divide 41 3 23)
    assertEqual "DIVIDE case 3" (7) (divide 1 3 10)

pointAddTest :: Assertion
pointAddTest = do
    assertEqual "ADD case 1" (x e) (x r)
    assertEqual "ADD case 2" (y e) (y r)
    where
    r = add(a, b, p, z)
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    z = 0
    a = Point 57812578399721466406643612476747524047439366610920293847719616669328532557630 26095420015187790947275343713856380659635149352631946998250298380930377583466
    b = Point 4032983015753143990395647783770666587927265353624430905763286836981504199392 44353125519324157186344456159742269880631179110473143840214086765587351124293
    e = Point 26935039991567691076591733863940144651417118571819271649648362736996102432831 53500484534674675258240380497240532959550103299561262531927294311330606411404

testTan :: Assertion
testTan = do
    assertEqual "TAN case 1" (r) (e)
    where
    r = tangent(a, z, p)
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    z = 0
    a = Point 87917229428110789366561422587307072970088695150214603900351294636804298290738 37579621872809717799779570009674914438024800886176540314162770583770981830863
    e = 101606632288299564871864232606624488455184839018758248336393796847926137794588

testDouble :: Assertion
testDouble = do
    assertEqual "DBL case 1" (x e) (x r)
    assertEqual "DBL case 2" (y e) (y r)
    where
    r = double(a, z, p)
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    z = 0
    a = Point 87917229428110789366561422587307072970088695150214603900351294636804298290738 37579621872809717799779570009674914438024800886176540314162770583770981830863
    e = Point 19277281477197177963613685635111727513957886411799201238917757645493897712993 847959926674921704613916930352312808004252888284294958523157455244708242291

testScalarMultiply :: Assertion
testScalarMultiply = do
   assertEqual "DBL case 1" (x e) (x d)
   assertEqual "DBL case 2" (y e) (y d)
    where
    d = pointMultiply(i, n, p, g, a)
    g = Point 55066263022277343669578718895168534326250603453777594175500187360389116729240 32670510020758816978083085130507043184471273380659243275938904335757337482424
    e = Point 23012785790775233531463825296277893910455674240231127288991410292022945764135 41599300880872949469217949279325578136298881846694700618797507821355387453499
    n = 42155799610086226451002635049826242127223948348387333047588118000654776762016
    p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    a = 0
    b = 7
    i = identify p

main :: IO ()
main = defaultMainWithOpts
       [
            (testCase "XGCD test" xgcdTest),
            (testCase "MOVINV test" modinvTest),
            (testCase "MUTPILY test" multiplyTest),
            (testCase "DIVIDE test" divideTest),
            (testCase "PADD test" pointAddTest),
            (testCase "TAN test" testTan),
            (testCase "DBL test" testDouble),
            (testCase "SCLR test" testScalarMultiply)
       ]
       mempty