#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ class std::vector<std::vector<int> >+;
#pragma link C++ class std::vector<std::vector<unsigned int> >+;
#pragma link C++ class std::vector<std::vector<short> >+;
#pragma link C++ class std::vector<std::vector<unsigned short> >+;
#pragma link C++ class std::vector<std::vector<long> >+;
#pragma link C++ class std::vector<std::vector<unsigned long> >+;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<double> >+;
#pragma link C++ class std::vector<std::string>+;
#pragma link C++ class std::vector<std::vector<std::string> >+;

#else
template class std::vector<std::vector<bool> >;
template class std::vector<std::vector<int> >;
template class std::vector<std::vector<double> >;
template class std::vector<std::vector<float> >;
#endif
