# femheart
   std::string dataname = "femheart.data";
   OptionsParser args(argc, argv);
   args.AddOption(&dataname, "-d", "--data", "Data file to use.");
   args.ParseCheck();
   
   std::vector<std::string> objectFilenames;
   //if (argc == 1)
   objectFilenames.push_back(dataname);
