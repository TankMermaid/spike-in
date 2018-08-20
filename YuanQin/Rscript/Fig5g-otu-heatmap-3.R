setwd("/mnt/bai/qinyuan/xiaoxuan/QA/")
print(paste("Your working directory is in",getwd()))

library(ggplot2)
library(plyr)
library(scales)
library(reshape2)
library(pheatmap)


################### AHWYJ #####################

data_AHWYJ <- read.csv("data/Fig5g-otu-heatmap-3-AHWYJ.csv", header = T)
all_AHWYJ <- gather(data_AHWYJ, key = group , value = logFC, `RAlogFC`:`AAlogFC`)

p_AHWYJ <- ggplot(all_AHWYJ, aes(group, AhWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "deepskyblue",mid="grey20",high = "yellow1", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("RAlogFC","AAlogFC"))+
  scale_y_discrete(limits=c("OTU_1029","OTU_104","OTU_108","OTU_111","OTU_1116","OTU_119","OTU_12","OTU_120","OTU_123","OTU_1279","OTU_137","OTU_143","OTU_146","OTU_17","OTU_176","OTU_18","OTU_180","OTU_185","OTU_189","OTU_229","OTU_23","OTU_234","OTU_237","OTU_259","OTU_262","OTU_269","OTU_279","OTU_281","OTU_283","OTU_289","OTU_291","OTU_294","OTU_2950","OTU_30","OTU_302","OTU_313","OTU_315","OTU_317","OTU_328","OTU_386","OTU_411","OTU_423","OTU_432","OTU_45","OTU_49","OTU_499","OTU_540","OTU_59","OTU_68","OTU_775","OTU_798","OTU_89","OTU_90","OTU_96","OTU_101","OTU_1027","OTU_1046","OTU_1095","OTU_115","OTU_1194","OTU_1230","OTU_130","OTU_132","OTU_1380","OTU_1390","OTU_140","OTU_145","OTU_154","OTU_166","OTU_167","OTU_179","OTU_181","OTU_187","OTU_190","OTU_20","OTU_206","OTU_209","OTU_21","OTU_2117","OTU_218","OTU_221","OTU_232","OTU_2486","OTU_261","OTU_270","OTU_278","OTU_285","OTU_292","OTU_296","OTU_298","OTU_304","OTU_308","OTU_31","OTU_311","OTU_312","OTU_319","OTU_321","OTU_331","OTU_342","OTU_363","OTU_370","OTU_371","OTU_374","OTU_38","OTU_387","OTU_40","OTU_410","OTU_416","OTU_418","OTU_43","OTU_443","OTU_445","OTU_447","OTU_450","OTU_486","OTU_498","OTU_506","OTU_513","OTU_515","OTU_54","OTU_558","OTU_564","OTU_565","OTU_630","OTU_679","OTU_703","OTU_745","OTU_757","OTU_765","OTU_767","OTU_80","OTU_83","OTU_9","OTU_92","OTU_931","OTU_937","OTU_948","OTU_99","OTU_611","OTU_1175","OTU_670","OTU_993","OTU_1093","OTU_122","OTU_124","OTU_134","OTU_14","OTU_1906","OTU_245","OTU_27","OTU_276","OTU_286","OTU_35","OTU_44","OTU_51","OTU_517","OTU_6","OTU_60","OTU_624","OTU_67","OTU_691","OTU_694","OTU_73","OTU_74","OTU_878","OTU_103","OTU_11","OTU_117","OTU_128","OTU_129","OTU_142","OTU_147","OTU_16","OTU_164","OTU_172","OTU_1729","OTU_184","OTU_205","OTU_208","OTU_22","OTU_226","OTU_236","OTU_244","OTU_248","OTU_266","OTU_290","OTU_3","OTU_307","OTU_309","OTU_326","OTU_346","OTU_36","OTU_403","OTU_412","OTU_414","OTU_42","OTU_47","OTU_520","OTU_53","OTU_621","OTU_701","OTU_71","OTU_720","OTU_811","OTU_84"))

p_AHWYJ

ggsave(paste("result/Fig5g-Heatmap-AHWYJ.pdf", sep=""), p_AHWYJ, width = 6, height = 8)

################### AHMH63 #####################

data_AHMH63 <- read.csv("data/Fig5g-otu-heatmap-3-AHMH63.csv", header = T)
all_AHMH63 <- gather(data_AHMH63, key = group , value = logFC, `RAlogFC`:`AAlogFC`)

p_AHMH63 <- ggplot(all_AHMH63, aes(group, AhMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "deepskyblue",mid="grey20",high = "yellow1", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("RAlogFC","AAlogFC"))+
  scale_y_discrete(limits=c("OTU_1380","OTU_143","OTU_185","OTU_214","OTU_23","OTU_237","OTU_259","OTU_267","OTU_30","OTU_302","OTU_371","OTU_374","OTU_443","OTU_499","OTU_68","OTU_92","OTU_96","OTU_104","OTU_1116","OTU_115","OTU_120","OTU_15","OTU_154","OTU_221","OTU_234","OTU_269","OTU_283","OTU_31","OTU_386","OTU_49","OTU_89","OTU_621","OTU_103","OTU_105","OTU_106","OTU_11","OTU_110","OTU_112","OTU_1134","OTU_114","OTU_116","OTU_117","OTU_1175","OTU_121","OTU_129","OTU_1310","OTU_133","OTU_140","OTU_1403","OTU_147","OTU_1565","OTU_16","OTU_165","OTU_167","OTU_170","OTU_171","OTU_173","OTU_174","OTU_175","OTU_178","OTU_183","OTU_187","OTU_188","OTU_19","OTU_1906","OTU_199","OTU_205","OTU_208","OTU_2117","OTU_215","OTU_219","OTU_22","OTU_220","OTU_224","OTU_227","OTU_230","OTU_231","OTU_2315","OTU_235","OTU_236","OTU_247","OTU_248","OTU_2588","OTU_27","OTU_270","OTU_276","OTU_278","OTU_279","OTU_286","OTU_290","OTU_292","OTU_299","OTU_3","OTU_307","OTU_308","OTU_324","OTU_331","OTU_354","OTU_357","OTU_363","OTU_38","OTU_384","OTU_39","OTU_393","OTU_400","OTU_408","OTU_412","OTU_414","OTU_43","OTU_444","OTU_454","OTU_47","OTU_470","OTU_486","OTU_491","OTU_50","OTU_506","OTU_51","OTU_517","OTU_520","OTU_526","OTU_53","OTU_543","OTU_558","OTU_564","OTU_57","OTU_6","OTU_60","OTU_63","OTU_635","OTU_64","OTU_654","OTU_665","OTU_67","OTU_694","OTU_707","OTU_71","OTU_736","OTU_751","OTU_791","OTU_84","OTU_91","OTU_936","OTU_94","OTU_107","OTU_124","OTU_125","OTU_128","OTU_138","OTU_160","OTU_172","OTU_184","OTU_186","OTU_211","OTU_226","OTU_244","OTU_245","OTU_266","OTU_297","OTU_301","OTU_309","OTU_326","OTU_35","OTU_42","OTU_44","OTU_643","OTU_686","OTU_720","OTU_73","OTU_811","OTU_82","OTU_878"))

p_AHMH63

ggsave(paste("result/Fig5g-Heatmap-AHMH63.pdf", sep=""), p_AHMH63, width = 6, height = 8)

################### HNWYJ #####################

data_HNWYJ <- read.csv("data/Fig5g-otu-heatmap-3-HNWYJ.csv", header = T)
all_HNWYJ <- gather(data_HNWYJ, key = group , value = logFC, `RAlogFC`:`AAlogFC`)

p_HNWYJ <- ggplot(all_HNWYJ, aes(group, HNWYJ)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "deepskyblue",mid="grey20",high = "yellow1", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("RAlogFC","AAlogFC"))+
  scale_y_discrete(limits=c("OTU_10","OTU_100","OTU_1095","OTU_110","OTU_115","OTU_119","OTU_120","OTU_123","OTU_131","OTU_1380","OTU_142","OTU_143","OTU_146","OTU_15","OTU_1523","OTU_154","OTU_1555","OTU_158","OTU_172","OTU_181","OTU_2","OTU_214","OTU_23","OTU_232","OTU_233","OTU_237","OTU_259","OTU_267","OTU_269","OTU_2795","OTU_283","OTU_290","OTU_294","OTU_2950","OTU_2952","OTU_30","OTU_302","OTU_31","OTU_312","OTU_33","OTU_34","OTU_36","OTU_374","OTU_411","OTU_412","OTU_443","OTU_444","OTU_492","OTU_498","OTU_499","OTU_54","OTU_61","OTU_62","OTU_66","OTU_68","OTU_72","OTU_80","OTU_89","OTU_9","OTU_90","OTU_92","OTU_98","OTU_99","OTU_1000","OTU_106","OTU_11","OTU_1107","OTU_1112","OTU_1121","OTU_1134","OTU_117","OTU_12","OTU_127","OTU_1310","OTU_1390","OTU_1403","OTU_144","OTU_151","OTU_1568","OTU_166","OTU_167","OTU_169","OTU_170","OTU_173","OTU_188","OTU_193","OTU_1991","OTU_206","OTU_2072","OTU_208","OTU_209","OTU_2117","OTU_212","OTU_218","OTU_226","OTU_229","OTU_247","OTU_256","OTU_268","OTU_2785","OTU_286","OTU_296","OTU_298","OTU_3","OTU_300","OTU_320","OTU_370","OTU_384","OTU_399","OTU_403","OTU_464","OTU_518","OTU_52","OTU_520","OTU_537","OTU_539","OTU_567","OTU_579","OTU_58","OTU_627","OTU_665","OTU_702","OTU_75","OTU_830","OTU_91","OTU_96","OTU_13","OTU_309","OTU_481","OTU_53","OTU_95","OTU_105","OTU_245","OTU_278","OTU_281","OTU_301","OTU_386","OTU_475","OTU_528","OTU_560","OTU_691","OTU_707","OTU_745","OTU_923","OTU_926","OTU_935","OTU_1029","OTU_103","OTU_1041","OTU_109","OTU_1230","OTU_128","OTU_129","OTU_133","OTU_134","OTU_137","OTU_138","OTU_16","OTU_160","OTU_161","OTU_162","OTU_165","OTU_182","OTU_1885","OTU_192","OTU_1944","OTU_1954","OTU_199","OTU_20","OTU_2077","OTU_217","OTU_22","OTU_225","OTU_24","OTU_248","OTU_2588","OTU_26","OTU_262","OTU_264","OTU_265","OTU_27","OTU_276","OTU_2775","OTU_28","OTU_285","OTU_289","OTU_311","OTU_343","OTU_346","OTU_35","OTU_354","OTU_356","OTU_357","OTU_38","OTU_405","OTU_416","OTU_423","OTU_438","OTU_45","OTU_46","OTU_49","OTU_502","OTU_515","OTU_541","OTU_57","OTU_59","OTU_592","OTU_638","OTU_67","OTU_675","OTU_680","OTU_694","OTU_77","OTU_775","OTU_78","OTU_826","OTU_83","OTU_931","OTU_101","OTU_102","OTU_107","OTU_118","OTU_122","OTU_132","OTU_14","OTU_141","OTU_153","OTU_17","OTU_18","OTU_19","OTU_21","OTU_210","OTU_213","OTU_222","OTU_231","OTU_266","OTU_29","OTU_37","OTU_373","OTU_4","OTU_41","OTU_51","OTU_55","OTU_574","OTU_6","OTU_60","OTU_64","OTU_65","OTU_670","OTU_69","OTU_7","OTU_81","OTU_87"))

p_HNWYJ

ggsave(paste("result/Fig5g-Heatmap-HNWYJ.pdf", sep=""), p_HNWYJ, width = 6, height = 8)


################### HNWYJ #####################

data_HNMH63 <- read.csv("data/Fig5g-otu-heatmap-3-HNMH63.csv", header = T)
all_HNMH63 <- gather(data_HNMH63, key = group , value = logFC, `RAlogFC`:`AAlogFC`)


p_HNMH63 <- ggplot(all_HNMH63, aes(group, HNMH63)) + 
  geom_tile(aes(fill = logFC)) + 
  scale_fill_gradient2(low = "deepskyblue",mid="grey20",high = "yellow1", midpoint = 0)+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))+
  scale_x_discrete(limits=c("RAlogFC","AAlogFC"))+
  scale_y_discrete(limits=c("OTU_10","OTU_100","OTU_1000","OTU_1095","OTU_115","OTU_119","OTU_123","OTU_131","OTU_1380","OTU_1403","OTU_142","OTU_143","OTU_146","OTU_15","OTU_150","OTU_154","OTU_1555","OTU_1568","OTU_166","OTU_1991","OTU_206","OTU_2072","OTU_214","OTU_23","OTU_233","OTU_237","OTU_2588","OTU_259","OTU_265","OTU_267","OTU_269","OTU_2795","OTU_283","OTU_290","OTU_294","OTU_2950","OTU_2952","OTU_3","OTU_30","OTU_302","OTU_31","OTU_312","OTU_33","OTU_36","OTU_374","OTU_411","OTU_492","OTU_498","OTU_499","OTU_518","OTU_564","OTU_62","OTU_72","OTU_80","OTU_89","OTU_9","OTU_90","OTU_92","OTU_96","OTU_98","OTU_99","OTU_11","OTU_1107","OTU_117","OTU_12","OTU_1310","OTU_173","OTU_181","OTU_188","OTU_208","OTU_213","OTU_217","OTU_226","OTU_231","OTU_245","OTU_256","OTU_264","OTU_2785","OTU_28","OTU_320","OTU_403","OTU_44","OTU_467","OTU_52","OTU_53","OTU_54","OTU_6","OTU_65","OTU_675","OTU_745","OTU_78","OTU_111","OTU_125","OTU_164","OTU_315","OTU_386","OTU_443","OTU_528","OTU_707","OTU_798","OTU_1029","OTU_106","OTU_1194","OTU_122","OTU_129","OTU_13","OTU_151","OTU_16","OTU_160","OTU_161","OTU_170","OTU_184","OTU_20","OTU_262","OTU_506","OTU_568","OTU_57","OTU_58","OTU_603","OTU_630","OTU_680","OTU_7","OTU_77","OTU_87","OTU_95","OTU_101","OTU_102","OTU_1047","OTU_110","OTU_1230","OTU_132","OTU_133","OTU_14","OTU_141","OTU_167","OTU_17","OTU_182","OTU_19","OTU_193","OTU_2","OTU_21","OTU_210","OTU_222","OTU_225","OTU_24","OTU_285","OTU_289","OTU_29","OTU_340","OTU_346","OTU_37","OTU_373","OTU_38","OTU_4","OTU_416","OTU_418","OTU_423","OTU_45","OTU_51","OTU_55","OTU_59","OTU_60","OTU_64","OTU_67","OTU_83","OTU_97"))

p_HNMH63

ggsave(paste("result/Fig5g-Heatmap-HNMH63.pdf", sep=""), p_HNMH63, width = 6, height = 8)











