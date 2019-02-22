READABLE = False

_data_ = {'MTSAT-1R': {'georange': (60., -60., 200., 80.), 'planck_const': {'ir': (9485.3, 1333.6, 0.34115, 0.99883), 'wv': (38792., 2132.7, 0.38020, 0.99911)},
                       'linear_const': {'ir':(-0.16552, 169.6753), 'wv':(-0.0316, 32.32898)}},
                  'MTSAT-2': {'georange': (60., -60., 205., 85.), 'planck_const': {'ir': (9474.4, 1333.1, 0.28349, 0.99903), 'wv': (38353., 2124.6, 0.36654, 0.99914)},
                              'linear_const':{'ir':(-0.16552, 169.6753), 'wv':(-0.0316, 32.32898)}},
                  'GOES-08': {'georange': (60., -60., 345., 225.), 'planck_const': {'ir': (9737.8, 1345.3, 0.37351, 0.99873), 'wv': (38792., 2132.7, 0.60603, 0.99858)},
                             'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-09': {'georange': (60., -60., 285., 165.), 'planck_const': {'ir': (9717.2, 1344.4, 0.36225, 0.99876), 'wv': (38732., 2131.6, 0.48409, 0.99887)},
                             'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-10': {'georange': (60., -60., 285., 165.), 'planck_const': {'ir': (9774.3, 1347.0, 0.27790, 0.99905), 'wv': (39086., 2138.1, 0.61437, 0.99857)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-11': {'georange': (60., -60., 285., 165.), 'planck_const': {'ir': (9653.4, 1341.5, 0.38284, 0.99869), 'wv': (38789., 2132.7, 0.59339, 0.99861)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-12': {'georange': (60., -60., 345., 225.), 'planck_const': {'ir': (9686.1, 1343.0, 0.37555, 0.99872), 'wv': (43703., 2219.2, 5.08328, 0.98872)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-13': {'georange': (60., -60., 345., 225.), 'planck_const': {'ir': (9805.0, 1348.5, 0.36348, 0.99876), 'wv': (42508., 2198.7, 3.96963, 0.99112)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-14': {'georange': (60., -60., 315., 195.), 'planck_const': {'ir': (9721.8, 1344.6, 0.36159, 0.99877), 'wv': (42172., 2192.9, 3.75614, 0.99154)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'GOES-15': {'georange': (60., -60., 285., 165.), 'planck_const': {'ir': (9765.6, 1346.6, 0.36151, 0.99877), 'wv': (42336., 2195.8, 3.75109, 0.99156)},
                              'linear_const':{'ir':(5.2285, 15.6854), 'wv':(38.8383, 29.1287)}},
                  'METEOSAT-8': {'alpha': {'band04': (0.9956), 'band05': (0.9962), 'band06': (0.9991), 'band07': (0.9996), 'band08': (0.9999), 'band09': (0.9983),
                                           'band10': (0.9988), 'band11': (0.9981)},
                                 'beta': {'band04': (3.410), 'band05': (2.218), 'band06': (0.478), 'band07': (0.179), 'band08': (0.060), 'band09': (0.625),
                                          'band10': (0.397), 'band11': (0.578)},
                                 'wnc': {'band04': (2567.330), 'band05': (1598.103), 'band06': (1362.081), 'band07': (1149.069), 'band08': (1034.343), 'band09': (930.647),
                                          'band10': (839.660), 'band11': (752.387)}},
                  'METEOSAT-9': {'alpha': {'band04': (0.9954), 'band05': (0.9963), 'band06': (0.9991), 'band07': (0.9996), 'band08': (0.9999), 'band09': (0.9983),
                                           'band10': (0.9988), 'band11': (0.9981)},
                                 'beta': {'band04': (3.438), 'band05': (2.185), 'band06': (0.470), 'band07': (0.179), 'band08': (0.056), 'band09': (0.640),
                                          'band10': (0.408), 'band11': (0.561)},
                                 'wnc': {'band04': (2568.832), 'band05': (1600.548), 'band06': (1360.330), 'band07': (1148.620), 'band08': (1035.289), 'band09': (931.700),
                                          'band10': (836.445), 'band11': (751.792)}},
                  'METEOSAT-10': {'alpha': {'band04': (0.9915), 'band05': (0.9960), 'band06': (0.9991), 'band07': (0.9996), 'band08': (0.9999), 'band09': (0.9983),
                                           'band10': (0.9988), 'band11': (0.9982)},
                                 'beta': {'band04': (2.9002), 'band05': (2.0337), 'band06': (0.4340), 'band07': (0.1714), 'band08': (0.0527), 'band09': (0.6084),
                                          'band10': (0.3882), 'band11': (0.5390)},
                                 'wnc': {'band04': (2547.771), 'band05': (1595.621), 'band06': (1360.377), 'band07': (1148.130), 'band08': (1034.715), 'band09': (929.842),
                                          'band10': (838.659), 'band11': (750.653)}},
                  'METEOSAT-11': {'alpha': {'band04': (0.9916), 'band05': (0.9959), 'band06': (0.9990), 'band07': (0.9996), 'band08': (0.9998), 'band09': (0.9983),
                                           'band10': (0.9988), 'band11': (0.9981)},
                                 'beta': {'band04': (2.9438), 'band05': (2.0780), 'band06': (0.4929), 'band07': (0.1731), 'band08': (0.0597), 'band09': (0.6256),
                                          'band10': (0.4002), 'band11': (0.5635)},
                                 'wnc': {'band04': (2555.280), 'band05': (1596.080), 'band06': (1361.748), 'band07': (1147.433), 'band08': (1034.851), 'band09': (931.122),
                                          'band10': (839.113), 'band11': (748.585)}}}       

#for segment calculation only
_hmw_data_for_seg_calc = {'EarthConst1':0.00669438, 'EarthConst2':0.99330562,
                          'EarthPolarRadius':6356.7523, 'SubLon':140.7, 'Distance':42164.,
                          'LOFF':2750.5, 'LFAC':20466275, 'COFF':2750.5, 'CFAC':20466275}
_FY4_data_for_seg_calc = {'EarthConst1':0.00669438, 'EarthConst2':0.99330562,
                          'EarthPolarRadius':6356.7523, 'SubLon':104.7, 'Distance':42164.,
                          'LOFF':2747.5, 'LFAC':20466274, 'COFF':2747.5, 'CFAC':20466274}  