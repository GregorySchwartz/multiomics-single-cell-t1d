#!/usr/bin/env stack
{- stack
  script
  --resolver lts-14.16
  --package system-filepath
  --package text
  --package turtle
  --package foldl
  --package async-pool
-}

{-# LANGUAGE OverloadedStrings #-}

import Turtle
import Data.Bool (bool)
import Data.Maybe (isJust, isNothing, fromJust)
import qualified Data.Text as T
import qualified Control.Foldl as Fold
import qualified Filesystem.Path.CurrentOS as FP
import qualified Control.Concurrent.Async.Pool as A

runCluster :: (Maybe Int, T.Text, T.Text, Maybe T.Text, T.Text, Maybe Int, Maybe Int, T.Text) -> IO ()
runCluster (cutSize, name, path, whitelist, norm, dropDim, rootCut, labelsFile) = sh $ do
    files <- reduce Fold.list . find (suffix ".csv") $ fromText path
    let fileArgs = concatMap (\x -> ["-m", format fp x]) files
        outputFile = "output_" <> name <> ".csv"
        outputCutFile Nothing Nothing = outputFile
        outputCutFile x rc = "output_" <> name <> maybe "" (\a -> "_" <> T.pack (show a) <> "_cut") x <> maybe "" (\a -> "_" <> T.pack (show a) <> "_rootCut") rc <> ".csv"
        outputCutFolder Nothing Nothing = "output_" <> name
        outputCutFolder x rc = "output_" <> name <> maybe "" (\a -> "_" <> T.pack (show a) <> "_cut") x <> maybe "" (\a -> "_" <> T.pack (show a) <> "_rootCut") rc

    liftIO $ print $ files

    output (FP.fromText $ outputCutFile cutSize rootCut)
        . inproc "too-many-cells" ( [ "make-tree"]
                                 <> bool fileArgs [] (isJust cutSize || isJust rootCut)
                                 <> [ "-o", outputCutFolder cutSize rootCut
                                    , "--draw-node-number"
                                    , "--draw-mark", "MarkModularity"
                                    , "--dendrogram-output", "tree_" <> (T.replace "--" "-" . T.replace "." "-" . T.replace "/" "-" $ labelsFile) <> ".pdf"
                                    , "--matrix-transpose"
                                    , "--draw-collection", "PieChart"
                                    , "--normalization", norm
                                    , "--filter-thresholds", "(1e-16, 1e-16)"
                                    ]
                                 <> bool [] [ "--prior", "output_" <> name] (isJust cutSize || isJust rootCut)
                                 <> maybe [] (\x -> [ "--smart-cutoff", T.pack . show $ x, "--min-size", "1" ]) cutSize
                                 <> bool [] [ "--labels-file", labelsFile ] (isJust cutSize || isJust rootCut)
                                 <> maybe [] (\x -> [ "--cell-whitelist-file", x ]) whitelist
                                 <> maybe ["--normalization", "TfIdfNorm"] (\x -> [ "--drop-dimension", "--lsa", T.pack $ show x ]) dropDim
                                 <> bool [] ([ "--root-cut", T.pack . show $ fromJust rootCut ]) (isJust rootCut && isNothing cutSize)
                                  )
        $ mempty

main :: IO ()
main = A.withTaskGroup 4 $ \g -> do
    _ <- A.mapConcurrently g runCluster $
        [
          -- Example
          (Just 5, "hpap_imc", "imc_data/", Nothing, "QuantileNorm", Nothing, Nothing, "labels_region_imc.csv")
        ]
    return ()
